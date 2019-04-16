import numpy
import skimulator.spline as spline
import scipy.sparse.linalg as linalg


class fitspline2d:
    def __init__(self, nx, ny, nxspline, nyspline, xcircle=True, ycircle=True):
        self.xcircle = xcircle
        self.ycircle = ycircle
        self.nx = nx
        self.ny = ny
        self.nxspline = nxspline
        self.nyspline = nyspline
        x = numpy.arange(nx) + 0.5
        y = numpy.arange(ny) + 0.5
        if self.xcircle == False:
            self.spx = spline.spline(3, nxspline - 2)
        else:
            self.spx = spline.spline_circular(3, nxspline)
        if self.ycircle == False:
            self.spy = spline.spline(3, nyspline - 2)
        else:
            self.spy = spline.spline_circular(3,nyspline)
        self.spx.init_uniform(0, nx)
        self.spy.init_uniform(0, ny)
        self.vx, self.ix, self.jx = self.spx.bspline_value(x)
        self.vy, self.iy, self.jy = self.spy.bspline_value(y)

        self.nspline=self.nxspline*self.nyspline

    def init_fit(self, doinvert=True):
        mat=numpy.zeros([self.nspline, self.nspline])

        for ii in range(self.nx):
            # print(ii)
            for jj in range(self.ny):
                for jjx in range(4):
                    for jjy in range(4):
                        for kkx in range(4):
                            for kky in range(4):
                                ixj = numpy.mod(self.ix[ii] + jjx, self.nxspline)
                                ixk = numpy.mod(self.ix[ii] + kkx, self.nxspline)
                                iyj = numpy.mod(self.iy[jj] + jjy, self.nyspline)
                                iyk = numpy.mod(self.iy[jj] + kky, self.nyspline)
                                ind0 = ixj + self.nxspline * iyj
                                ind1 = ixk + self.nxspline * iyk
                                mat[ind0, ind1] += (self.vx[kkx,ii]
                                                    * self.vy[kky,jj]
                                                    * self.vx[jjx,ii]
                                                    * self.vy[jjy,jj])
        if self.xcircle == False:
            for ii in range(self.nyspline):
                mat[ii * self.nxspline, ii * self.nxspline] += 1
            for ii in range(self.nyspline):
                mat[self.nxspline-1+ii*self.nxspline, ii*self.nxspline] += 1
        if self.ycircle == False:
            for ii in range(self.nxspline):
                mat[ii, ii] += 1
            for ii in range(self.nxspline):
                ind0 = ii+(self.nyspline-1)*self.nxspline
                ind1 = (self.nyspline-1)*self.nxspline
                mat[ind0, ind1] += 1
        if doinvert is True:
            self.imat=numpy.linalg.pinv(mat)
            self.inv = True
        else:
            self.imat = mat
            self.inv = False

    def fit(self,data):

        nspline=self.nspline
        vec=numpy.zeros([nspline])

        for ii in range(self.nx):
            for jj in range(self.ny):
                for jjx in range(4):
                    for jjy in range(4):
                        ixj = numpy.mod(self.ix[ii] + jjx, self.nxspline)
                        iyj = numpy.mod(self.iy[jj] + jjy, self.nyspline)
                        vec[ixj + self.nxspline * iyj] += self.vx[jjx, ii]*self.vy[jjy, jj] * data[ii, jj]

        if self.xcircle == False:
            for ii in range(self.nyspline):
                vec[ii*self.nxspline]+=data[0,ii]
            for ii in range(self.nyspline):
                vec[self.nxspline-1+ii*self.nxspline]+=data[self.nx-1,ii]
        if self.ycircle == False:
            for ii in range(self.nxspline):
                vec[ii]+=data[ii,0]
            for ii in range(self.nxspline):
                vec[ii+(self.nyspline-1)*self.nxspline]+=data[ii,self.ny-1]
        if self.inv == False:
            self.wres,infoinv=linalg.cg(self.imat,vec)
        else:
            self.wres=numpy.dot(self.imat, vec)

    def getparam(self):
        return self.wres

    def setparam(self,wres):
        self.wres=wres

    def transform(self, x, y):
        n = x.shape[0]
        res = numpy.zeros([n])

        vx, ix, jx = self.spx.bspline_value(x)
        vy, iy, jy = self.spy.bspline_value(y)

        self.wres=self.wres.reshape(self.nyspline, self.nxspline)

        for jj in range(4):
            for ii in range(4):
                iyj = numpy.mod(iy + jj, self.nyspline)
                ixi = numpy.mod(ix + ii, self.nxspline)
                res[:] += vy[jj, :] * vx[ii, :] * self.wres[iyj, ixi]

        return res


class ted_tas():
    def __init__(self, wres, nxspline, nyspline):
        self.mysp2={}
        self.n_time = 1200
        self.n_az = 120
        for i in range(3):
            self.mysp2[i]=fitspline2d(self.n_time, self.n_az, nxspline, nyspline)
            self.mysp2[i].setparam(wres[0,:])


    def transform(self,t_orbit,az,t_year):
        '''
        #====================================================================
        # t_orbit : orbital time inside [0,1]; 0 begining of the orbit, 1 end
        of the orbit (if t_orbit <0 or t_orbit>1 then BUG!
        # az  : azimuth in degree
        # t_year : date in year inside [0,1]; 0 begining of the year,
                                              1 end of the year
        '''
        # Number of points for the time dimension
        n_time = self.n_time
        az = numpy.mod(az, 360)
        # Discretisation step for azimuth
        delta_az = 360 / self.n_az
        fit_1 = self.mysp2[0].transform(t_orbit * n_time, az / delta_az)
        cos_year = numpy.cos(t_year * 2 * numpy.pi)
        fit_t = self.mysp2[1].transform(t_orbit * n_time, az / delta_az)
        fit_2 = fit_t * numpy.cos(t_year * 2 * numpy.pi)
        fit_3 = fit_t * numpy.sin(t_year * 2 * numpy.pi)
        return (fit_1 + fit_2 + fit_3)
