cimport numpy as cnp
import numpy as np

import thermodyn

'''
Code signature for the emanuel convection fortran code
'''

cdef extern:
    void emanuel_convection(double *temp, double *q, double *qs,
                            double *u, double *v, double *tracers,
                            double *pmid, double *pint, int *nlevs,
                            int *max_conv_lev, int *num_tracers,
                            double *dt, int *conv_state,
                            double *dtemp, double *dq, double *du,
                            double *dv, double *dtracers, double *precip,
                            double *downdraft_vel_scale, double *downdraft_temp_scale,
                            double *downdraft_q_scale, double *cloud_base_mass_flux,
                            double *cape)


'''
Arrays and variables to pass to the fortran code
'''

cdef cnp.double_t[::1] temp, q, qs, u, v, dtemp, dq, du, dv, pmid, pint
cdef cnp.double_t[::1,:] tracers, dtracers, old_cloud_base_mass_flux, cape_conv
cdef cnp.int_t[::1,:] convection_state

cdef double dt, precip, downdraft_vel_scale
cdef double downdraft_temp_scale, downdraft_q_scale, cloud_base_mass_flux, cape

cdef int nlevs, num_tracers, max_conv_lev, nlats, nlons, initialised


def init_emanuel_convection(num_levs, num_lats, num_lons, time_step, num_trac=0):

    global old_cloud_base_mass_flux, nlevs, num_tracers, tracers, dtracers,\
        convection_state, dt, nlats, nlons, initialised, max_conv_lev,\
        temp, q, qs, u, v, pmid, pint, dtemp, du, dv, dq, cape_conv

    initialised = 0
    # Two arrays needed by the emanuel scheme:
    # old_cloud_base_mass_flux holds the mass flux from the previous call
    # convection_state hold the status of convection in each column

    cape_conv = np.zeros((num_lons, num_lats), dtype=np.double, order='F')
    old_cloud_base_mass_flux = np.zeros((num_lons, num_lats), dtype=np.double, order='F')
    convection_state = np.zeros((num_lons, num_lats), dtype=np.int, order='F')

    if num_trac == 0:
        num_tracers = 1
    else:
        num_tracers = num_trac

    dt = time_step

    # Copy grid dimensions
    nlevs = num_levs
    nlons = num_lons
    nlats = num_lats
    max_conv_lev = nlevs - 3 #Max value = nlevs - 1


    #Initialise arrays
    temp = np.zeros(nlevs, dtype=np.double, order='F')
    q = np.zeros(nlevs, dtype=np.double, order='F')
    qs = np.zeros(nlevs, dtype=np.double, order='F')
    u = np.zeros(nlevs, dtype=np.double, order='F')
    v = np.zeros(nlevs, dtype=np.double, order='F')
    pmid = np.zeros(nlevs, dtype=np.double, order='F')
    pint = np.zeros(nlevs+1, dtype=np.double, order='F')
    dtemp = np.zeros(nlevs, dtype=np.double, order='F')
    du = np.zeros(nlevs, dtype=np.double, order='F')
    dv = np.zeros(nlevs, dtype=np.double, order='F')
    dq = np.zeros(nlevs, dtype=np.double, order='F')


    dtracers = np.zeros((num_tracers, nlevs), dtype=np.double, order='F')
    tracers = np.zeros((num_tracers, nlevs), dtype=np.double, order='F')

    print 'Finished initialising emanuel convection scheme'
    
    initialised = 1


def driver(cnp.ndarray[cnp.double_t, ndim=3] t_ext,
           cnp.ndarray[cnp.double_t, ndim=3] q_ext,
           cnp.ndarray[cnp.double_t, ndim=3] u_ext,
           cnp.ndarray[cnp.double_t, ndim=3] v_ext,
           cnp.ndarray[cnp.double_t, ndim=3] p_ext,
           cnp.ndarray[cnp.double_t, ndim=3] pint_ext,
           cnp.ndarray[cnp.double_t, ndim=4] tra_ext = None):

    cdef double temp_cbmf
    cdef int conv_state
    cdef  cnp.double_t[::1] ptr

    global old_cloud_base_mass_flux, nlevs, num_tracers, tracers, dtracers,\
        convection_state, dt, nlats, nlons, initialised, cloud_base_mass_flux,\
        temp, q, qs, u, v, pmid, pint, dtemp, du, dv, dq, max_conv_lev,\
        precip, downdraft_q_scale, downdraft_temp_scale, downdraft_vel_scale, cape, cape_conv

    if initialised == 0:
        raise ValueError, 'Emanuel scheme not initialised.'

    arr_shape = np.array(t_ext).shape
    #create output arrays
    temp_tend = np.empty(arr_shape, dtype=np.double, order='F')
    q_tend = np.empty(arr_shape, dtype=np.double, order='F')
    u_tend = np.empty(arr_shape, dtype=np.double, order='F')
    v_tend = np.empty(arr_shape, dtype=np.double, order='F')
    precipitation = np.empty((nlons,nlats), dtype=np.double, order='F')

    sat_q = np.array(thermodyn.qs(np.asfortranarray(t_ext), np.asfortranarray(p_ext)/100.)/1000., order='C')

    # Convert to Numpy Type for multiplication with memoryviews
    time_step = np.double(dt)

    t_ext = np.asarray(t_ext, order='C')
    q_ext = np.asarray(q_ext, order='C')
    u_ext = np.asarray(u_ext, order='C')
    v_ext = np.asarray(v_ext, order='C')
    p_ext = np.asarray(p_ext, order='C')
    pint_ext = np.asarray(pint_ext, order='C')

    for lon in range(nlons):
        for lat in range(nlats):

            #Initialise all arrays for convect
            cloud_base_mass_flux = old_cloud_base_mass_flux[lon,lat]
            cape = 0.

            ptr = t_ext[lon,lat,:]
            temp[:] = ptr[:]

            ptr =  sat_q[lon,lat,:]
            qs[:] = ptr[:]

            ptr = q_ext[lon,lat,:]
            q[:] = ptr[:]

            ptr = u_ext[lon,lat,:]
            u[:] = ptr[:]

            ptr = v_ext[lon,lat,:]
            v[:] = ptr[:]

            #Scheme expects pressures in hPa
            ptr = p_ext[lon,lat,:]/100.
            pmid[:] = ptr[:]

            ptr = pint_ext[lon,lat,:]/100.
            pint[:] = ptr[:]

            #Call fortran code
            emanuel_convection(<double *>&temp[0], <double *>&q[0],
                               <double *>&qs[0], <double *>&u[0],
                               <double *>&v[0], <double *>&tracers[0,0],
                               <double *>&pmid[0], <double *>&pint[0],

                               &nlevs, &max_conv_lev, &num_tracers,
                               &dt, &conv_state,

                               <double *>&dtemp[0], <double *>&dq[0],
                               <double *>&du[0], <double *>&dv[0],
                               <double *>&dtracers[0,0],

                               &precip, &downdraft_vel_scale,
                               &downdraft_temp_scale,
                               &downdraft_q_scale, &cloud_base_mass_flux, &cape)

            # Copy to output arrays
            temp_tend[lon,lat,:] = np.asfortranarray(dtemp[:])*dt
            u_tend[lon,lat,:] = np.asfortranarray(du[:])*dt
            v_tend[lon,lat,:] = np.asfortranarray(dv[:])*dt
            q_tend[lon,lat,:] = np.asfortranarray(dq[:])*dt
            precipitation[lon,lat] = precip

            old_cloud_base_mass_flux[lon,lat] = cloud_base_mass_flux
            convection_state[lon,lat] = conv_state
            cape_conv[lon,lat] = cape

            Tdot = np.empty(arr_shape, dtype=np.double, order='F')
            Tdot[lon,lat,:] = np.asfortranarray(dtemp[:])

            Qdot = np.empty(arr_shape, dtype=np.double, order='F')
            Qdot[lon,lat,:] = np.asfortranarray(dq[:])
    
    cape_array = np.asfortranarray(cape_conv)
    #return [temp_tend, q_tend, u_tend, v_tend, precipitation]
    return [temp_tend, q_tend, precipitation, Tdot, Qdot, cape_array]


