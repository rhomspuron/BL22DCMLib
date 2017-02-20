import math

EV2REVA = 0.2624682843
c=299792458
h=4.13566727e-15
aSi=5.43102088e-10 # Simo macro su spec
hc = 12398.419 #eV *Angstroms
PMAC_OVERFLOW = 8388608 #2**23 encoder register 24 bits overflow



cdef double energy2bragg (double energy, double vcm_pitch, double d,
                   double offset):
    """
    Method to calculate the bragg angle
    :param energy: Energy in eV
    :param vcm_pitch_rad: VCM pitch angle in rad
    :param d: constant d of the crystal
    :param offset: offset of the crystal
    :return: bragg angle in degree
    """
    cdef double bragg_rad, bragg_deg

    try:
        bragg_rad = math.asin(hc/2/d/energy) + 2 * vcm_pitch + offset
    except ZeroDivisionError as e:
        bragg_rad = float('nan')

    bragg_deg = math.degrees(bragg_rad)
    return bragg_deg


cdef double bragg2energy (double bragg, double vcm_pitch, double d,
                   double offset):
    """
    Method to calculate the bragg angle
    :param bragg: Bragg angle in degree
    :param vcm_pitch_rad: VCM pitch angle in rad
    :param d: constant d of the crystal
    :param offset: offset of the crystal
    :return: Energy in eV
    """

    cdef double bragg_rad, energy

    bragg_rad = math.radians(bragg)
    try:
        energy = hc / (2 * d * math.sin(bragg_rad - 2 * vcm_pitch - offset))
    except ZeroDivisionError as e:
        energy = float('nan')

    return energy



cdef long bragg2encoder(double bragg, double bragg_spu, double bragg_offset, 
                         long pmac_offset, long pmac_enc):
    """
    Function to calculate the encoder value for one bragg angle.
    :param bragg: Bragg angle in degree
    :param bragg_spu: Bragg step per unit (Sardana)
    :param bragg_offset: Bragg offset (Sardana)
    :param pmac_offset: Pmac internal offset of the bragg motor
    :param pmac_enc: Pmac internal value of the encoder
    :return:
    """


    cdef:
        long braggMotorOffsetEncCounts, enc_value
        double offset

    braggMotorOffsetEncCounts = <long> (bragg_offset * bragg_spu)
    offset = pmac_offset - pmac_enc + braggMotorOffsetEncCounts

    enc_value = <long> ((bragg * bragg_spu) - offset)

    if enc_value > PMAC_OVERFLOW:
        enc_value = enc_value - 2 * PMAC_OVERFLOW
    elif enc_value < -PMAC_OVERFLOW:
        enc_value = enc_value + 2 * PMAC_OVERFLOW

    return enc_value


cdef double encoder2bragg(long encoder, double bragg_spu, double bragg_offset,
                          long pmac_offset, long pmac_enc):
    """
    Function to calcule the encoder value for one bragg angle.
    """


    #translations from raw counts to degrees
    #getting an offset between position and encoder register (offset = 2683367)

    cdef:
        long braggMotorOffsetEncCounts, enc_value
        double offset

    braggMotorOffsetEncCounts = <long> (bragg_offset * bragg_spu)
    offset = pmac_offset - pmac_enc + braggMotorOffsetEncCounts

    enc_value = encoder + offset
    if enc_value > PMAC_OVERFLOW:
        enc_value = enc_value - 2 * PMAC_OVERFLOW
    elif enc_value < -PMAC_OVERFLOW:
        enc_value = enc_value + 2 * PMAC_OVERFLOW

    return enc_value/bragg_spu


cdef long energy2encoder (double energy, double vcm_pitch_rad, double d,
                          double offset, double bragg_spu,
                          double bragg_offset, long pmac_offset, long pmac_enc):

    cdef:
        double bragg
        long encoder

    bragg = energy2bragg(energy, vcm_pitch_rad, d, offset)
    encoder = bragg2encoder(bragg, bragg_spu, bragg_offset, pmac_offset,
                            pmac_enc)
    return encoder


cdef double encoder2energy (long encoder, double vcm_pitch_rad, double d,
                          double offset, double bragg_spu,
                          double bragg_offset, long pmac_offset, long pmac_enc):

    cdef double bragg, energy
    bragg = encoder2bragg(encoder,bragg_spu, bragg_offset, pmac_offset,
                          pmac_enc)
    energy = bragg2energy(bragg, vcm_pitch_rad, d, offset)
    return energy



cdef enegies4encoders(long* encoders, double vcm_pitch_rad, double d,
                      double offset, double bragg_spu,
                      double bragg_offset, long pmac_offset, long pmac_enc):

    pass
