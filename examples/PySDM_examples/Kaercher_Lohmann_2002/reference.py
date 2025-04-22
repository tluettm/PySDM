import numpy as np
from PySDM.physics.constants import si

def critical_supersaturation(temperature):
    return 2.349 - temperature / 259.

# crit thres
def critical_supersaturation_spich2023 (t):
    s20=1.67469
    s21=0.00228125
    s22=-1.36989e-05
    return s20+s21*t+s22*t*t


def bulk_model_reference(initial_temperature, updraft=0.1):
    n_hom_ice = None

    if initial_temperature == 220.:
        if updraft == 0.1:
            n_hom_ice = 148121.413358197
        if updraft == 1.:
            n_hom_ice = 7268664.77542974

    if initial_temperature == 216.:
        if updraft == 1.:
            n_hom_ice = 10475282.894692907



    return( n_hom_ice / si.metre ** 3  )

def bulk_model_reference_array():

    dsd_list = np.array([0,1,2])
    initial_temperatures = np.array([196., 216., 236.])
    updrafts = np.array([0.05, 0.1, 0.3, 0.5, 1., 3., 5., 10.])

    dim_size = ( np.shape(dsd_list)[0], np.shape(initial_temperatures)[0], np.shape(updrafts)[0] )
    ni_bulk_ref = np.zeros(dim_size)


    # DSD 1 T = 196
    ni_bulk_ref[0,0,0] = 643686.1316903427
    ni_bulk_ref[0,0,1] = 2368481.0609527444
    ni_bulk_ref[0,0,2] = 20160966.984670535
    ni_bulk_ref[0,0,3] = 49475281.81718969
    ni_bulk_ref[0,0,4] = 131080662.23620115
    ni_bulk_ref[0,0,5] = 401046528.70428866
    ni_bulk_ref[0,0,6] = 627442148.3402529
    ni_bulk_ref[0,0,7] = 1151707310.2210448

    # DSD 1 T = 216
    ni_bulk_ref[0,1,0] = 60955.84292640147
    ni_bulk_ref[0,1,1] = 189002.0792186534
    ni_bulk_ref[0,1,2] = 1200751.6897658105
    ni_bulk_ref[0,1,3] = 2942110.815055958
    ni_bulk_ref[0,1,4] = 10475282.894692907
    ni_bulk_ref[0,1,5] = 90871045.40856971
    ni_bulk_ref[0,1,6] = 252175505.460412
    ni_bulk_ref[0,1,7] = 860335156.4717773

    # DSD 1 T = 236
    ni_bulk_ref[0,2,0] = 13049.108886452004
    ni_bulk_ref[0,2,1] = 40422.244759544985
    ni_bulk_ref[0,2,2] = 237862.49854786208
    ni_bulk_ref[0,2,3] = 545315.7805748513
    ni_bulk_ref[0,2,4] = 1707801.469906006
    ni_bulk_ref[0,2,5] = 11128055.66932415
    ni_bulk_ref[0,2,6] = 27739585.111447476
    ni_bulk_ref[0,2,7] = 101799566.47225031

    ########################

    # DSD 2 T = 196
    ni_bulk_ref[1,0,0] = 637413.
    ni_bulk_ref[1,0,1] = 2344493.
    ni_bulk_ref[1,0,2] = 19966877.
    ni_bulk_ref[1,0,3] = 49069527.
    ni_bulk_ref[1,0,4] = 130283347.
    ni_bulk_ref[1,0,5] = 399266774.
    ni_bulk_ref[1,0,6] = 624750911.
    ni_bulk_ref[1,0,7] = 1146655515.

    # DSD 2 T = 216
    ni_bulk_ref[1,1,0] = 60403.
    ni_bulk_ref[1,1,1] = 187283.
    ni_bulk_ref[1,1,2] = 1189550.
    ni_bulk_ref[1,1,3] = 2914046.
    ni_bulk_ref[1,1,4] = 10371219.
    ni_bulk_ref[1,1,5] = 89923341.
    ni_bulk_ref[1,1,6] = 249676049.
    ni_bulk_ref[1,1,7] = 853573871.

    # DSD 2 T = 236
    ni_bulk_ref[1,2,0] = 12912.
    ni_bulk_ref[1,2,1] = 40018.
    ni_bulk_ref[1,2,2] = 235542.
    ni_bulk_ref[1,2,3] = 540005.
    ni_bulk_ref[1,2,4] = 1691068.
    ni_bulk_ref[1,2,5] = 11015461.
    ni_bulk_ref[1,2,6] = 27451520.
    ni_bulk_ref[1,2,7] = 100693311.

    ########################

    # DSD 3 T = 196
    ni_bulk_ref[2,0,0] = 656510.
    ni_bulk_ref[2,0,1] = 2417540.
    ni_bulk_ref[2,0,2] = 20556685.
    ni_bulk_ref[2,0,3] = 50299515.
    ni_bulk_ref[2,0,4] = 132694835.
    ni_bulk_ref[2,0,5] = 404648346.
    ni_bulk_ref[2,0,6] = 632892416.
    ni_bulk_ref[2,0,7] = 1161943841.

    # DSD 3 T = 216
    ni_bulk_ref[2,1,0] = 62083
    ni_bulk_ref[2,1,1] = 192509
    ni_bulk_ref[2,1,2] = 1223618
    ni_bulk_ref[2,1,3] = 2999419
    ni_bulk_ref[2,1,4] = 10687903
    ni_bulk_ref[2,1,5] = 92807057
    ni_bulk_ref[2,1,6] = 257269740
    ni_bulk_ref[2,1,7] = 874045873

    # DSD 3 T = 236
    ni_bulk_ref[2,2,0] = 13339.
    ni_bulk_ref[2,2,1] = 41258.
    ni_bulk_ref[2,2,2] = 242613.
    ni_bulk_ref[2,2,3] = 556177.
    ni_bulk_ref[2,2,4] = 1742003.
    ni_bulk_ref[2,2,5] = 11358229.
    ni_bulk_ref[2,2,6] = 28328686.
    ni_bulk_ref[2,2,7] = 104063397.


    return ni_bulk_ref