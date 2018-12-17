import os
import zlib
import base64
import struct
import hashlib
import logging

logger = logging.getLogger(__name__)


def _text_file_md5(file_path, block_size=1024*1024):
    """"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(file_path)

    md5 = hashlib.md5()
    with open(file_path, 'r', encoding='utf-8') as f:
        while True:
            data = f.read(block_size)
            if not data:
                break
            md5.update(data.encode('utf-8'))
    file_hash =  md5.hexdigest()
    return file_hash


def get_hash_from_grid_params(p):
    """"""
    orbit_path = p.filesat  # Path to a text file containing the elevation, the cycle and the track of the satellite
    orbit_cols = p.order_orbit_col  # None|int[3] (order of columns in orbit file)
    modelbox = p.modelbox  # None|float[4] (lonmin, lonmax, latmin, latmax)
    rotation_speed = p.rotation_speed  # float (rotation speed of the antenna in tr/min)
    cycle = p.cycle  # float (beam cycle duration)
    list_pos = p.list_pos  # float[N] (position of each beam)
    list_angle = p.list_angle  # float[N] (inclination of each beam)
    list_shift = p.list_shift  # int[N] (time shift wrt nadir for each beam)
    shift_lon = p.shift_lon  # None|float (longitude shift in degrees)
    shift_time = p.shift_time  # None|float (time shift in days)

    orbit_file_hash = _text_file_md5(orbit_path)
    orbit_md5 = orbit_file_hash.encode('utf-8')
    logger.debug('Orbit file MD5: {}'.format(orbit_file_hash))

    # Flags for optional parameters
    has_cols = orbit_cols is not None
    has_modelbox = modelbox is not None
    has_shift_lon = shift_lon is not None
    has_shift_time = shift_time is not None

    fmt = '!B32sffI{}{}{}'.format('f' * len(list_pos),
                                  'f' * len(list_angle),
                                  'i' * len(list_shift))
    logger.debug('Mandatory params serialization format: {}'.format(fmt))

    flags = 0b00000000
    opts = []
    if has_cols:
        opts.append(struct.pack('!BBB', orbit_cols))
    else:
        flags = flags + 0b00000001

    if has_modelbox:
        opts.append(struct.pack('!ffff', modelbox))
    else:
        flags = flags + 0b00000010

    if has_shift_lon:
        opts.append(struct.pack('!f', shift_lon))
    else:
        flags = flags + 0b00000100

    if has_shift_time:
        opts.append(struct.pack('!f', shift_time))
    else:
        flags = flags + 0b00001000

    h = struct.pack(fmt, flags, orbit_md5, rotation_speed, cycle,
                    len(list_pos), *list_pos, *list_angle, *list_shift)
    for d in opts:
        h = h + d

    return h


def get_grid_params_from_hash(h):
    """"""
    flags, orbit_md5, rot, cycle, nbeams = struct.unpack('!B32sffI', h[:45])
    has_cols = 0 >= (flags & 0b00000001)
    has_modelbox = 0 >= (flags & 0b00000010)
    has_shift_lon = 0 >= (flags & 0b00000100)
    has_shift_time = 0 >= (flags & 0b00001000)

    istart = 45
    iend = 45 + 4 * nbeams  # assuming 1 float <=> 4 bytes
    list_pos = struct.unpack('!{}f'.format(nbeams), h[istart:iend])
    istart = iend

    iend = istart + 4 * nbeams  # assuming 1 float <=> 4 bytes
    list_angle = struct.unpack('!{}f'.format(nbeams), h[istart:iend])
    istart = iend

    iend = istart + 4 * nbeams  # assuming 1 int <=> 4 bytes
    list_shift = struct.unpack('!{}i'.format(nbeams), h[istart:iend])
    istart = iend

    if has_cols:
        iend = istart + 3
        orbit_cols = struct.unpack('!BBB', h[istart:iend])
        istart = iend
    else:
        orbit_cols = None

    if has_modelbox:
        iend = istart + 4 * 4 # assuming 1 float <=> 4 bytes
        modelbox = struct.unpack('!ffff', h[istart:iend])
        istart = iend
    else:
        modelbox = None

    if has_shift_lon:
        iend = istart + 4  # assuming 1 float <=> 4 bytes
        shift_lon = struct.unpack('!f', h[istart:iend])[0]
        istart = iend
    else:
        shift_lon = None

    if has_shift_time:
        iend = istart + 4  # assuming 1 float <=> 4 bytes
        shift_time = struct.unpack('!f', h[istart:iend])[0]
        istart = iend
    else:
        shift_time = None

    result = {'order_orbit_col': orbit_cols,
              'rotation_speed': rot,
              'cycle': cycle,
              'list_pos': list_pos,
              'list_angle': list_angle,
              'list_shift': list_shift,
              'modelbox': modelbox,
              'shift_lon': shift_lon,
              'shift_time': shift_time,
              'orbit_file_md5': orbit_md5}

    logger.debug('')
    logger.debug('=== Reconstructed parameters ===')
    logger.debug('orbit file MD5: {}'.format(orbit_md5.decode('utf-8')))
    logger.debug('rotation_speed: {}'.format(rot))
    logger.debug('cycle: {}'.format(cycle))
    logger.debug('nbeams: {}'.format(nbeams))
    logger.debug('list_pos: {}'.format(list_pos))
    logger.debug('list_angle: {}'.format(list_angle))
    logger.debug('list_shift: {}'.format(list_shift))
    logger.debug('order_orbit_col: {}'.format(orbit_cols))
    logger.debug('modelbox: {}'.format(modelbox))
    logger.debug('shift_lon: {}'.format(shift_lon))
    logger.debug('shift_time: {}'.format(shift_time))

    return result


def get_b64_gzipped_hash(p):
    """"""
    h = get_hash_from_grid_params(p)
    compressed = zlib.compress(h, 9)
    b64_bytes = base64.b64encode(compressed)
    b64_str = b64_bytes.decode('utf-8')
    return b64_str


def revert_b64_gzipped_hash(b64_str):
    """"""
    b64_bytes = b64_str.encode('utf-8')
    compressed = base64.b64decode(b64_bytes)
    h = zlib.decompress(compressed)
    p = get_grid_params_from_hash(h)
    return p


if '__main__' == __name__:
    import sys
    import logging
    import params_8beams_fram as p

    main_logger = logging.getLogger()
    main_logger.handlers = []
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    main_logger.addHandler(handler)
    main_logger.setLevel(logging.DEBUG)

    try:
        h = compute_grid_hash(p)
    except FileNotFoundError:
        _, e, _ = sys.exc_info()
        logger.error('Orbit file "{}" not found'.format(e.args[0]))

    # Compress information as much as possible...
    print(base64.b64encode(zlib.compress(h, 9)).decode('utf-8'))

    # Maybe use the CRC32 in the filename and include the hash in the NetCDF
    # attributes
    print(zlib.crc32(h))

    read_grid_params_from_hash(h)
