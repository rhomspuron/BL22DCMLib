import math

EV2REVA = 0.2624682843
c=299792458
h=4.13566727*10**(-15)
aSi=5.43102088*10**(-10) # Simo macro su spec
hc = 12398.419 #eV *Angstroms
PMAC_OVERFLOW = 8388608 #2**23 encoder register 24 bits overflow



cdef double energy2bragg(double energy, str crystal, double vcm_pitch_mrad, 
                         double dSi111=3.1354161, double offsetSi111=0, 
                         double dSi311=1.637418, double offsetSi311=0):
     """
     Function to calculate the DCM bragg angle for one energy.
     """
     cdef double vcm_pitch_rad, d, offset, bragg_rad

     vcm_pitch_rad = vcm_pitch_mrad / 1000
     if crystal.lower() == 'si111':
         d = dSi111
         offset = offsetSi111
     elif crystal.lower() == 'si311':
         d = dSi311
         offset = offsetSi311
        
     try:
         bragg_rad = math.asin(hc/2/d/energy) + 2 * vcm_pitch_rad + offset
     except ZeroDivisionError,e:
         bragg_rad = float('nan')
     bragg_deg = math.degrees(bragg_rad)
     return bragg_deg

cdef long bragg2encoder(double bragg, double bragg_spu, double bragg_offset, 
                         long pmac_offset, long pmac_enc):
        """
        Function to calcule the encoder value for one bragg angle.
        """

        # Translations from  degrees to raw counts. Getting an offset between
        # position and encoder register (offset = 2683367)
       
        cdef long braggMotorOffsetEncCounts, enc_value
        cdef double offset

        braggMotorOffsetEncCounts = <long> (bragg_offset * bragg_spu)
        offset = pmac_offset - pmac_enc + braggMotorOffsetEncCounts

        enc_value = <long> ((bragg * bragg_spu) - offset)
        
        if enc_value > PMAC_OVERFLOW:
            enc_value = enc_value - 2 * PMAC_OVERFLOW
        elif enc_value < -PMAC_OVERFLOW:
            enc_value = enc_value + 2 * PMAC_OVERFLOW
        
        return enc_value
        

cpdef get_enc_table(double start_energy, double end_energy, 
                    long nr_points, double int_time, double bragg_spu, 
                    double bragg_offset, long pmac_offset, 
                    long pmac_enc, str crystal, 
                    double vcm_pitch_mrad, double dSi111=3.1354161, 
                    double offsetSi111=0, double dSi311=1.637418, 
                    double offsetSi311=0):

    cdef double th1, th2, scan_time, v_theta
    th1= energy2bragg(start_energy, crystal, vcm_pitch_mrad, dSi111, 
                      offsetSi111, dSi311, offsetSi311)
    th2= energy2bragg(end_energy, crystal, vcm_pitch_mrad, dSi111, 
                      offsetSi111, dSi311, offsetSi311)
   
    scan_time = nr_points * int_time
        
    v_theta = (th1-th2)/scan_time
    enc_table = []
    for i in range(nr_points):
        theta = (-v_theta*int_time*i+ th1)
        enc = bragg2encoder(theta, bragg_spu, bragg_offset, pmac_offset, 
                            pmac_enc,)
        enc_table.append(enc)
