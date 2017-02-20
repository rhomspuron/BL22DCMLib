import math

EV2REVA = 0.2624682843
c=299792458
h=4.13566727*10**(-15)
aSi=5.43102088*10**(-10) # Simo macro su spec
hc = 12398.419 #eV *Angstroms
PMAC_OVERFLOW = 8388608 #2**23 encoder register 24 bits overflow



def energy2bragg(energy, crystal, vcm_pitch_mrad, dSi111=3.1354161, 
                  offsetSi111=0, dSi311=1.637418, offsetSi311=0):
     """
     Function to calculate the DCM bragg angle for one energy.
     """
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

def bragg2encoder(bragg, bragg_spu, bragg_offset, pmac_offset, pmac_enc):
        """
        Function to calcule the encoder value for one bragg angle.
        """

        # Translations from  degrees to raw counts. Getting an offset between
        # position and encoder register (offset = 2683367)
       
        braggMotorOffsetEncCounts = bragg_offset * bragg_spu
        offset = pmac_offset - pmac_enc + braggMotorOffsetEncCounts

        enc_value = (bragg * bragg_spu) - offset
        
        if enc_value > PMAC_OVERFLOW:
            enc_value = enc_value - 2 * PMAC_OVERFLOW
        elif enc_value < -PMAC_OVERFLOW:
            enc_value = enc_value + 2 * PMAC_OVERFLOW
        
        return enc_value
        

def get_enc_table(start_energy, end_energy, nr_points, int_time, bragg_spu, 
                  bragg_offset, pmac_offset, pmac_enc, crystal, vcm_pitch_mrad, 
                  dSi111=3.1354161, offsetSi111=0, dSi311=1.637418, offsetSi311=0 ):
    
    th1= energy2bragg(start_energy, crystal, vcm_pitch_mrad, dSi111, 
                      offsetSi111, dSi311, offsetSi311)
    th2= energy2bragg(end_energy,crystal, vcm_pitch_mrad, dSi111, 
                      offsetSi111, dSi311, offsetSi311)
   
    scan_time = nr_points * int_time
        
    v_theta = (th1-th2)/scan_time
    enc_table = []
    for i in range(nr_points):
        theta = (-v_theta*int_time*i+ th1)
        enc = bragg2encoder(theta, bragg_spu, bragg_offset, pmac_offset, pmac_enc)
        enc_table.append(enc)



 
