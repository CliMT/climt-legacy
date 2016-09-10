cimport numpy as cnp
import numpy as np

# Simple physics package fortran code signature
# For a description, see Reed and Jablonowski (2012)
cdef extern:
   void dcmipSimplePhysics(int *numColumns, int *numLevels,
                           double *dt, double *latitude,
                           double *t, double *q,
                           double *u, double* v,
                           double *pressure, double *ifacePressure,
                           double *layerThickness, double *invThickness,
                           double *surfPressure, double *precipRate,
                           int *test)


#Arrays to pass to fortran. The dimensions are (levels, num_columns) for 3d and (num_columns) for 2d
cdef cnp.double_t[::1] temp, q, u, v, p, p_iface, thickness, inv_thickness
cdef double surf_press, precip, latitude
cdef int num_cols, num_levs, initialised, cyclone_simulation
cdef double time_step

def init_arrays(int columns, int lons, int lats, int levels, double dt, int simulate_cyclone):
    '''
    Initialise global arrays and some constants
    '''

    global temp, q, u, v, p, p_iface, thickness, inv_thickness,\
        surf_press, precip, latitude, num_cols, num_levs, initialised, time_step

    temp = np.empty((levels), dtype=np.double, order='F')
    q = np.empty((levels), dtype=np.double, order='F')
    u = np.empty((levels), dtype=np.double, order='F')
    v = np.empty((levels), dtype=np.double, order='F')
    p = np.empty((levels), dtype=np.double, order='F')
    thickness = np.empty((levels), dtype=np.double, order='F')
    inv_thickness = np.empty((levels), dtype=np.double, order='F')

    p_iface = np.empty((levels+1), dtype=np.double, order='F')

    test = np.zeros((lons,lats,levels), dtype=np.double, order='F')

    #surf_press = np.empty(columns, dtype=np.double, order='F')
    #precip = np.empty(columns, dtype=np.double, order='F')
    #latitude = np.empty(columns, dtype=np.double, order='F')

    #TODO Until we are happy and stable, only one column
    num_cols = 1
    num_levs = levels
    time_step = dt
    cyclone_simulation = simulate_cyclone

    initialised = 1

#Returns a list of tendencies
def get_tendencies(cnp.ndarray[cnp.double_t, ndim=3] u_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] v_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] temp_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] p_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] p_iface_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] q_ext,
                    cnp.ndarray[cnp.double_t, ndim=2] surf_press_ext,
                    cnp.ndarray[cnp.double_t, ndim=1] latitude_ext):
    '''
    This function takes fetches tendencies from the Simple Physics code. It
    

    '''

    global temp, q, u, v, p, p_iface, latitude, thickness, inv_thickness,\
        surf_press, precip, num_cols, num_levs, initialised, time_step

    cdef  cnp.double_t[::1] ptr
    if not initialised: 
        print 'init_arrays not called!'
        return

    nlons, nlats = np.asfortranarray(surf_press_ext).shape
    u_out = np.zeros((nlons,nlats,num_levs),order='F')
    v_out = np.zeros((nlons,nlats,num_levs),order='F')
    t_out = np.zeros((nlons,nlats,num_levs),order='F')
    q_out = np.zeros((nlons,nlats,num_levs),order='F')
    precip_out = np.zeros((nlons,nlats),order='F')
    #thickness = np.empty(temp_ext.shape, dtype=np.double, order='F')
    #inv_thickness = np.empty(temp_ext.shape, dtype=np.double, order='F')
    for lon in range(nlons):
        for lat in range(nlats):

            ptr = np.asfortranarray(u_ext[lon,lat,::-1])
            u[:] = ptr[:]

            ptr = np.asfortranarray((p_iface_ext[lon,lat,0:-1] - p_iface_ext[lon,lat,1::])[::-1])
            thickness[:] = ptr[:]
            
            ptr = np.asfortranarray(temp_ext[lon,lat,::-1])
            temp[:] = ptr[:]

            ptr = np.asfortranarray(v_ext[lon,lat,::-1])
            v[:] = ptr[:]

            ptr = np.asfortranarray(q_ext[lon,lat,::-1])
            q[:] = ptr[:]

            ptr = np.asfortranarray(p_ext[lon,lat,::-1])
            p[:] = ptr[:]

            ptr = np.asfortranarray(p_iface_ext[lon,lat,::-1])
            p_iface[:] = ptr[:]

            ptr = 1./np.asfortranarray(thickness)
            inv_thickness[:] = ptr[:] 

            surf_press = surf_press_ext[lon,lat]

            latitude = latitude_ext[lat]


            #Call fortran code with these arguments
            dcmipSimplePhysics(&num_cols, &num_levs,
                       &time_step, &latitude,
                       <double *>&temp[0], <double *>&q[0],
                       <double *>&u[0], <double *>&v[0],
                       <double *>&p[0], <double *>&p_iface[0],
                       <double *>&thickness[0], <double *>&inv_thickness[0],
                       &surf_press, &precip,
                       &cyclone_simulation)

            t_out[lon,lat,:] = (temp[::-1] - temp_ext[lon,lat,:])
            u_out[lon,lat,:] = (u[::-1] - u_ext[lon,lat,:])
            v_out[lon,lat,:] = (v[::-1] - v_ext[lon,lat,:])
            q_out[lon,lat,:] = (q[::-1] - q_ext[lon,lat,:])
            precip_out[lon,lat] = precip
    
    '''
    print 'T has nan: ', np.any(np.isnan(t_out))
    print 'U has nan: ', np.any(np.isnan(u_out))
    print 'V has nan: ', np.any(np.isnan(v_out))
    print 'Q has nan: ', np.any(np.isnan(q_out))
    print 'P has nan: ', np.any(np.isnan(precip_out))

    print 'T tendency: ', (np.max(t_out))
    print 'U tendency: ', (np.max(u_out))
    print 'V tendency: ', (np.max(v_out))
    print 'Q tendency: ', (np.max(q_out))
    print 'P tendency: ', (np.max(precip_out))
    '''

    return t_out, u_out, v_out, q_out, precip_out

