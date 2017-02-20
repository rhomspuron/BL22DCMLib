import math
from cython cimport boundscheck, wraparound
from cython.parallel cimport prange
import numpy as np


EV2REVA = 0.2624682843
c=299792458
h=4.13566727e-15
aSi=5.43102088e-10 # Simo macro su spec
hc = 12398.419 #eV *Angstroms
PMAC_OVERFLOW = 8388608 #2**23 encoder register 24 bits overflow



cpdef double energy2bragg (double energy, double vcm_pitch, double d,
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


cpdef double bragg2energy (double bragg, double vcm_pitch, double d,
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



cpdef long bragg2encoder(double bragg, double bragg_spu, double bragg_offset,
                         long pmac_offset, long pmac_enc):
    """
    Function to calculate the encoder value for one bragg angle.
    :param bragg: Bragg angle in degrees
    :param bragg_spu: Bragg step per unit (Sardana)
    :param bragg_offset: Bragg offset (Sardana)
    :param pmac_offset: Pmac internal offset of the bragg motor
    :param pmac_enc: Pmac internal value of the encoder
    :return: encoder value
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


cpdef double encoder2bragg(long encoder, double bragg_spu, double bragg_offset,
                          long pmac_offset, long pmac_enc):
    """
    Function to calcule the bragg angle value for one encoder.

    :param encoder: Encoder value
    :param bragg_spu: Bragg step per unit (Sardana)
    :param bragg_offset: Bragg offset (Sardana)
    :param pmac_offset: Pmac internal offset of the bragg motor
    :param pmac_enc: Pmac internal value of the encoder
    :return: Bragg angle in degrees
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


cpdef long energy2encoder (double energy, double vcm_pitch_rad, double d,
                          double offset, double bragg_spu,
                          double bragg_offset, long pmac_offset, long pmac_enc):
    """
    Function to convert from energy to encoder value

    :param energy: Energy in eV
    :param vcm_pitch_rad: VCM pitch angle in rad
    :param d: constant d of the crystal
    :param offset: offset of the crystal
    :param bragg_spu: Bragg step per unit (Sardana)
    :param bragg_offset: Bragg offset (Sardana)
    :param pmac_offset: Pmac internal offset of the bragg motor
    :param pmac_enc: Pmac internal value of the encoder
    :return: encoder value
    """

    cdef:
        double bragg
        long encoder

    bragg = energy2bragg(energy, vcm_pitch_rad, d, offset)
    encoder = bragg2encoder(bragg, bragg_spu, bragg_offset, pmac_offset,
                            pmac_enc)
    return encoder


cpdef double encoder2energy (long encoder, double vcm_pitch_rad, double d,
                          double offset, double bragg_spu,
                          double bragg_offset, long pmac_offset, long pmac_enc):

    """
    Function to convert from encoder to energy value

    :param encoder: encoder value
    :param vcm_pitch_rad: VCM pitch angle in rad
    :param d: constant d of the crystal
    :param offset: offset of the crystal
    :param bragg_spu: Bragg step per unit (Sardana)
    :param bragg_offset: Bragg offset (Sardana)
    :param pmac_offset: Pmac internal offset of the bragg motor
    :param pmac_enc: Pmac internal value of the encoder
    :return: Energy in eV
    """

    cdef double bragg, energy
    bragg = encoder2bragg(encoder,bragg_spu, bragg_offset, pmac_offset,
                          pmac_enc)
    energy = bragg2energy(bragg, vcm_pitch_rad, d, offset)
    return energy



@boundscheck(False)
@wraparound(False)
cpdef enegies4encoders(long[:] encoders, double vcm_pitch_rad, double d,
                       double offset, double bragg_spu,
                       double bragg_offset, long pmac_offset, long pmac_enc):

    """
    Funcion to transform an encoder array to energies array

    :param encoders: numpy array of encoder value
    :param vcm_pitch_rad: VCM pitch angle in rad
    :param d: constant d of the crystal
    :param offset: offset of the crystal
    :param bragg_spu: Bragg step per unit (Sardana)
    :param bragg_offset: Bragg offset (Sardana)
    :param pmac_offset: Pmac internal offset of the bragg motor
    :param pmac_enc: Pmac internal value of the encoder
    :return: numpy array with energies values in eV
    """
    cdef:
        long i, N
        double[:] values

    N = encoders.shape[0]
    values = np.zeros(N, dtype=np.float64)

    for i in prange(N, nogil=True):
        values[i] = encoder2energy(encoders[i], vcm_pitch_rad, d, offset,
                                   bragg_spu, bragg_offset, pmac_offset,
                                   pmac_enc)
    return np.asarray(values)

@boundscheck(False)
@wraparound(False)
cpdef enegies4encoders(double[:] energies, double vcm_pitch_rad, double d,
                       double offset, double bragg_spu,
                       double bragg_offset, long pmac_offset, long pmac_enc):

    """
    Funcion to transform an encoder array to energies array

    :param encoders: numpy array of encoder value
    :param vcm_pitch_rad: VCM pitch angle in rad
    :param d: constant d of the crystal
    :param offset: offset of the crystal
    :param bragg_spu: Bragg step per unit (Sardana)
    :param bragg_offset: Bragg offset (Sardana)
    :param pmac_offset: Pmac internal offset of the bragg motor
    :param pmac_enc: Pmac internal value of the encoder
    :return: numpy array with energies values in eV
    """
    cdef:
        long i, N
        double[:] values

    N = energies.shape[0]
    values = np.zeros(N, dtype=np.int32)

    for i in prange(N, nogil=True):
        values[i] = energy2encoder(energies[i], vcm_pitch_rad, d, offset,
                                   bragg_spu, bragg_offset, pmac_offset,
                                   pmac_enc)
    return np.asarray(values)

