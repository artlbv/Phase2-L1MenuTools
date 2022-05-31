from array import *

off = {}
offrate = {}
onl = {}
onlrate = {}


off['tkMuon'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0, 60.0, 63.0, 66.0, 69.0, 72.0, 75.0])
offrate['tkMuon'] =  array('d', [11230.9, 9820.3, 1230.6, 272.4, 94.2, 42.4, 23.6, 13.9, 8.4, 5.1, 3.9, 3.2, 2.1, 1.8, 1.5, 1.3, 1.1, 1.0, 0.9, 0.9, 0.8, 0.8, 0.8, 0.6, 0.5, 0.5])

onl['tkMuon'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0, 60.0, 63.0, 66.0, 69.0, 72.0, 75.0])
onlrate['tkMuon'] =  array('d', [11230.9, 5012.5, 719.7, 182.4, 65.6, 34.2, 19.0, 10.5, 6.7, 4.3, 3.3, 2.3, 1.9, 1.5, 1.4, 1.1, 1.0, 0.9, 0.9, 0.8, 0.8, 0.8, 0.5, 0.5, 0.4, 0.4])

off['tkMuonBarrel'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
offrate['tkMuonBarrel'] =  array('d', [677.1, 661.1, 314.8, 91.6, 32.9, 15.6, 9.7, 6.4, 3.8, 2.1, 1.6, 1.4, 1.0, 0.8, 0.8, 0.8, 0.6, 0.6, 0.6, 0.6])

onl['tkMuonBarrel'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
onlrate['tkMuonBarrel'] =  array('d', [677.1, 538.6, 221.1, 65.0, 23.2, 13.0, 8.2, 4.7, 2.8, 1.8, 1.4, 1.1, 0.8, 0.8, 0.8, 0.6, 0.6, 0.6, 0.6, 0.6])

off['tkMuonOverlap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
offrate['tkMuonOverlap'] =  array('d', [1137.6, 1099.2, 250.3, 48.3, 15.7, 7.7, 4.5, 2.5, 1.8, 1.1, 0.9, 0.8, 0.6, 0.6, 0.5, 0.3, 0.3, 0.1, 0.1, 0.1])

onl['tkMuonOverlap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
onlrate['tkMuonOverlap'] =  array('d', [1137.6, 843.8, 139.2, 30.3, 10.7, 6.5, 3.3, 2.1, 1.5, 0.9, 0.8, 0.6, 0.6, 0.5, 0.4, 0.3, 0.1, 0.1, 0.1, 0.1])

off['tkMuonEndcap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
offrate['tkMuonEndcap'] =  array('d', [10094.5, 8637.5, 685.7, 134.2, 45.7, 19.1, 9.3, 5.0, 2.8, 1.9, 1.4, 1.0, 0.6, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3])

onl['tkMuonEndcap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
onlrate['tkMuonEndcap'] =  array('d', [10094.5, 3859.3, 366.7, 88.0, 31.7, 14.6, 7.5, 3.7, 2.4, 1.6, 1.0, 0.6, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1])

off['standaloneMuon'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0, 60.0, 63.0, 66.0, 69.0, 72.0, 75.0])
offrate['standaloneMuon'] =  array('d', [14991.8, 13574.9, 2648.7, 1089.7, 445.8, 190.8, 86.3, 56.5, 42.1, 33.4, 29.3, 21.2, 19.0, 15.1, 14.0, 13.4, 11.3, 11.0, 9.8, 9.3, 7.9, 7.4, 6.4, 6.0, 5.9, 5.8])

onl['standaloneMuon'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0, 60.0, 63.0, 66.0, 69.0, 72.0, 75.0])
onlrate['standaloneMuon'] =  array('d', [14991.8, 5101.2, 1156.7, 361.7, 160.2, 73.9, 50.5, 35.6, 30.5, 21.8, 19.3, 15.2, 12.5, 11.9, 10.4, 10.0, 8.6, 7.0, 6.5, 6.2, 6.1, 5.4, 5.2, 5.1, 4.3, 4.3])

off['standaloneMuonBarrel'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
offrate['standaloneMuonBarrel'] =  array('d', [1281.9, 1281.9, 1271.3, 618.0, 220.9, 72.6, 28.7, 16.2, 12.2, 8.4, 6.4, 4.8, 4.1, 3.5, 3.2, 2.9, 2.6, 2.6, 2.4, 2.1])

onl['standaloneMuonBarrel'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
onlrate['standaloneMuonBarrel'] =  array('d', [1281.9, 1281.9, 414.9, 123.3, 46.0, 20.3, 14.1, 9.9, 7.0, 5.4, 4.4, 3.7, 3.3, 3.0, 2.8, 2.5, 2.5, 2.1, 2.1, 2.0])

off['standaloneMuonOverlap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
offrate['standaloneMuonOverlap'] =  array('d', [780.0, 780.0, 464.5, 191.0, 109.3, 62.1, 28.0, 21.2, 17.0, 14.0, 13.9, 9.4, 9.4, 7.3, 7.2, 7.2, 5.9, 5.9, 5.0, 5.0])

onl['standaloneMuonOverlap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
onlrate['standaloneMuonOverlap'] =  array('d', [780.0, 777.8, 313.2, 98.8, 57.4, 26.9, 20.6, 14.2, 14.0, 9.4, 9.4, 7.3, 5.9, 5.9, 5.0, 5.0, 4.1, 3.2, 3.2, 3.2])

off['standaloneMuonEndcap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
offrate['standaloneMuonEndcap'] =  array('d', [13936.3, 12419.1, 1032.7, 305.3, 125.1, 60.6, 31.7, 20.8, 14.4, 11.7, 9.6, 7.1, 5.7, 4.4, 3.7, 3.4, 2.8, 2.6, 2.5, 2.2])

onl['standaloneMuonEndcap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
onlrate['standaloneMuonEndcap'] =  array('d', [13936.3, 3359.9, 467.9, 148.1, 60.6, 28.8, 17.6, 12.3, 10.2, 7.2, 5.7, 4.3, 3.4, 3.0, 2.6, 2.5, 2.1, 1.8, 1.3, 1.1])

off['tkMuonStub'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0, 60.0, 63.0, 66.0, 69.0, 72.0, 75.0])
offrate['tkMuonStub'] =  array('d', [28416.9, 21838.3, 1338.7, 300.9, 107.2, 49.1, 26.6, 16.2, 10.3, 6.9, 5.2, 4.3, 3.1, 2.8, 2.3, 2.1, 1.7, 1.4, 1.4, 1.3, 1.3, 1.3, 1.3, 0.9, 0.9, 0.9])

onl['tkMuonStub'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0, 60.0, 63.0, 66.0, 69.0, 72.0, 75.0])
onlrate['tkMuonStub'] =  array('d', [28416.9, 9979.5, 806.4, 205.2, 74.4, 38.9, 21.5, 12.5, 8.5, 5.5, 4.4, 3.2, 2.8, 2.3, 2.2, 1.7, 1.4, 1.4, 1.3, 1.3, 1.3, 1.3, 0.9, 0.9, 0.9, 0.9])

off['tkMuonStubBarrel'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
offrate['tkMuonStubBarrel'] =  array('d', [530.9, 530.9, 282.8, 97.4, 36.7, 17.1, 10.4, 6.9, 4.3, 2.7, 2.4, 2.0, 1.4, 1.3, 1.1, 1.0, 0.8, 0.8, 0.8, 0.8])

onl['tkMuonStubBarrel'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
onlrate['tkMuonStubBarrel'] =  array('d', [530.9, 530.9, 205.1, 67.4, 25.3, 13.9, 8.8, 5.0, 3.5, 2.5, 2.1, 1.4, 1.3, 1.1, 1.1, 0.8, 0.8, 0.8, 0.8, 0.8])

off['tkMuonStubOverlap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
offrate['tkMuonStubOverlap'] =  array('d', [2705.2, 1972.1, 204.2, 43.8, 14.9, 7.2, 4.6, 2.5, 1.9, 1.1, 0.8, 0.7, 0.5, 0.5, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1])

onl['tkMuonStubOverlap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
onlrate['tkMuonStubOverlap'] =  array('d', [2705.2, 708.3, 114.5, 27.5, 10.0, 6.0, 3.3, 2.1, 1.5, 0.9, 0.8, 0.6, 0.5, 0.4, 0.4, 0.3, 0.1, 0.1, 0.1, 0.1])

off['tkMuonStubEndcap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
offrate['tkMuonStubEndcap'] =  array('d', [28187.9, 21167.1, 869.7, 161.7, 55.9, 24.9, 11.7, 6.8, 4.1, 3.1, 2.0, 1.6, 1.2, 1.0, 0.8, 0.8, 0.6, 0.6, 0.6, 0.5])

onl['tkMuonStubEndcap'] =  array('d', [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0, 57.0])
onlrate['tkMuonStubEndcap'] =  array('d', [28187.9, 9162.9, 494.1, 111.3, 39.3, 19.1, 9.4, 5.4, 3.5, 2.2, 1.6, 1.2, 1.0, 0.8, 0.8, 0.6, 0.6, 0.6, 0.5, 0.5])


##EG rates

off['tkElectron'] =  array('d', [10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0, 70.0, 73.0, 76.0, 79.0, 82.0, 85.0, 88.0, 91.0, 94.0, 97.0])
offrate['tkElectron'] =  array('d', [884.4, 570.6, 334.0, 201.8, 112.9, 76.9, 55.2, 40.5, 27.3, 22.5, 17.4, 14.0, 11.8, 9.5, 7.4, 6.2, 5.3, 4.2, 3.5, 3.2, 2.4, 2.2, 2.1, 1.8, 1.7, 1.6, 1.4, 1.2, 1.1, 1.0])

off['tkElectronBarrel'] = array('d', [10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0])
offrate['tkElectronBarrel'] = array('d', [620.1, 364.9, 210.0, 134.5, 74.4, 51.1, 38.9, 28.8, 19.3, 15.9, 12.2, 10.1, 8.4, 6.9, 5.2, 4.5, 3.8, 3.5, 2.9, 2.6])

off['tkElectronEndcap'] = array('d', [10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0])
offrate['tkElectronEndcap'] = array('d', [271.6, 210.2, 125.8, 68.3, 38.7, 25.9, 16.5, 11.7, 7.9, 6.7, 5.3, 3.9, 3.4, 2.6, 2.1, 1.7, 1.4, 0.7, 0.6, 0.6])

off['tkIsoElectron'] =  array('d', [10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0, 70.0, 73.0, 76.0, 79.0, 82.0, 85.0, 88.0, 91.0, 94.0, 97.0])
offrate['tkIsoElectron'] =  array('d', [470.6, 281.0, 151.7, 92.9, 51.8, 35.7, 27.2, 20.5, 14.5, 12.3, 10.6, 8.6, 7.5, 6.4, 5.5, 4.7, 3.6, 3.0, 2.3, 2.0, 1.6, 1.4, 1.4, 1.3, 1.1, 1.1, 0.9, 0.8, 0.6, 0.5])

off['standaloneElectron'] =  array('d', [10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0, 70.0, 73.0, 76.0, 79.0, 82.0, 85.0, 88.0, 91.0, 94.0, 97.0])
offrate['standaloneElectron'] =  array('d', [21719.7, 19089.8, 5527.7, 1862.3, 887.7, 431.0, 264.0, 173.0, 114.4, 78.7, 58.4, 44.8, 34.6, 28.0, 22.7, 18.2, 14.0, 11.3, 8.7, 7.0, 5.4, 4.7, 4.2, 3.6, 3.3, 3.1, 2.7, 2.3, 1.9, 1.4])

off['tkPhotonIso'] =  array('d', [10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0, 70.0, 73.0, 76.0, 79.0, 82.0, 85.0, 88.0, 91.0, 94.0, 97.0])
offrate['tkPhotonIso'] =  array('d', [8825.4, 3482.5, 978.1, 474.4, 243.0, 150.3, 104.5, 76.1, 52.0, 39.4, 31.0, 24.8, 20.5, 17.4, 13.4, 10.7, 8.3, 6.3, 5.1, 3.8, 3.2, 3.0, 2.7, 2.6, 2.4, 2.0, 1.6, 1.3, 0.9, 0.8])

off['tkPhotonIsoBarrel'] =  array('d', [10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0, 70.0, 73.0, 76.0, 79.0, 82.0, 85.0, 88.0, 91.0, 94.0, 97.0])
offrate['tkPhotonIsoBarrel'] =  array('d', [1868.9, 825.1, 463.4, 296.4, 158.6, 103.5, 74.4, 55.8, 37.3, 27.9, 23.0, 18.3, 15.3, 13.1, 10.1, 8.6, 6.9, 5.2, 4.2, 3.3, 2.8, 2.6, 2.3, 2.3, 2.0, 1.6, 1.4, 1.0, 0.6, 0.6])

off['tkPhotonIsoEndcap'] =  array('d', [10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0, 70.0, 73.0, 76.0, 79.0, 82.0, 85.0, 88.0, 91.0, 94.0, 97.0])
offrate['tkPhotonIsoEndcap'] =  array('d', [7400.1, 2736.2, 524.2, 180.7, 85.5, 47.2, 30.4, 20.5, 14.8, 11.7, 8.2, 6.7, 5.4, 4.5, 3.4, 2.3, 1.5, 1.1, 0.9, 0.6, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3])

###hadronic ttbar
#off['puppiPhase1HT'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0])
#offrate['puppiPhase1HT'] =  array('d', [31038.0, 5536.7, 4481.5, 1790.8, 1099.2, 629.6, 372.8, 226.0, 143.9, 94.2, 66.4, 47.5, 32.7, 23.5, 17.4, 13.8, 10.5, 8.5, 6.9, 4.8, 4.0, 3.3, 2.6, 2.0, 1.8, 1.4, 1.3, 1.1, 0.8, 0.6, 0.6, 0.6, 0.6, 0.5, 0.4, 0.4, 0.3, 0.3])


off['puppiMET'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0])
offrate['puppiMET'] =  array('d', [31038.0, 13791.6, 2921.4, 620.2, 156.9, 48.3, 18.1, 9.2, 5.3, 3.4, 2.2, 1.8, 1.5, 1.3, 1.0, 0.7, 0.6, 0.5])

off['puppiPhase1Jet'] =  array('d', [40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0])
offrate['puppiPhase1Jet'] =  array('d', [16324.4, 5363.9, 2085.5, 847.5, 393.2, 205.5, 115.6, 70.1, 45.7, 28.8, 18.3, 12.5, 8.2, 6.2, 4.9, 3.6, 2.6, 2.1, 1.8, 1.5])

off['puppiPhase1JetExt'] =  array('d', [40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0])
offrate['puppiPhase1JetExt'] =  array('d', [17492.9, 6541.9, 2936.1, 1401.1, 670.8, 345.8, 198.3, 119.9, 78.9, 51.1, 33.5, 23.9, 17.5, 13.8, 11.1, 8.9, 7.2, 6.2, 5.5, 5.0])

#taueta 2.4
off['NNPuppiTauLoose'] =  array('d', [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0])
offrate['NNPuppiTauLoose'] =  array('d', [2406.3, 2182.7, 1846.7, 1449.3, 1049.6, 726.3, 525.7, 387.3, 294.1, 224.5, 179.0, 141.8, 115.2, 93.5, 79.0, 69.1, 58.7, 50.9, 44.7, 39.4, 35.3, 31.0, 28.4, 25.8, 23.4, 21.3, 19.1, 17.6, 16.3, 15.4])

#2.1 #these aresomehow faulty?
#off['NNPuppiTauLoose'] =  array('d', [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0])
#offrate['NNPuppiTauLoose'] =  array('d', [2162.6, 1950.4, 1638.8, 1278.8, 924.7, 644.1, 471.0, 349.4, 266.4, 204.1, 163.3, 130.0, 106.4, 87.0, 73.1, 64.7, 54.9, 47.2, 41.5, 36.6, 32.9, 29.2, 26.8, 24.4, 22.4, 20.4, 18.3, 17.1, 15.8, 14.9])

off['NNPuppiTauLoose21'] =  array('d', [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0])
offrate['NNPuppiTauLoose21'] =  array('d', [2095.7, 1888.7, 1583.8, 1234.8, 891.5, 621.9, 456.3, 339.7, 259.8, 199.5, 160.1, 127.2, 104.1, 84.8, 71.2, 62.9, 53.3, 45.8, 40.1, 35.2, 31.9, 28.2, 25.9, 23.5, 21.7, 19.8, 17.8, 16.5, 15.2, 14.4])

off['NNPuppiTauLooseBarrel'] =  array('d', [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0])
offrate['NNPuppiTauLooseBarrel'] =  array('d', [1599.6, 1413.7, 1152.2, 878.1, 623.3, 451.1, 342.8, 259.7, 202.3, 156.9, 127.7, 101.4, 84.5, 69.7, 59.0, 52.7, 44.3, 38.4, 33.8, 29.9, 26.9, 23.9, 22.0, 20.3, 18.9, 17.3, 15.5, 14.3, 13.3, 12.6])

#hadronic ttbar
#off['trackerHT'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0])
#offrate['trackerHT'] =  array('d', [27857.3, 17890.9, 9170.4, 4480.2, 2252.9, 1212.2, 693.2, 420.3, 266.3, 179.7, 123.8, 88.4, 64.4, 48.1, 36.7, 28.2, 22.2, 18.4, 15.1, 12.4, 10.3, 8.2, 7.4, 5.9, 5.0, 4.2, 3.3, 2.8, 2.5, 1.8, 1.4, 1.3, 1.2, 1.1, 1.1, 0.9, 0.9, 0.7])
#onl['trackerHT'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0])
#onlrate['trackerHT'] =  array('d', [2809.5, 586.8, 182.1, 73.1, 34.2, 18.8, 11.3, 7.0, 4.4, 2.6, 1.4, 1.1, 0.9, 0.6, 0.6, 0.5, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1])

off['caloHT'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0])
offrate['caloHT'] =  array('d', [20275.1, 15173.6, 12163.9, 9764.8, 7278.0, 5558.7, 4265.6, 3188.0, 2398.7, 1801.8, 1339.2, 991.8, 726.1, 529.6, 383.7, 278.1, 200.6, 141.6, 99.9, 68.6, 49.0, 34.2, 23.1, 16.4, 10.6, 8.1, 5.5, 4.0, 2.9, 2.3, 1.6, 1.3, 0.9, 0.8, 0.7, 0.6, 0.6, 0.5])
onl['caloHT'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0])
onlrate['caloHT'] =  array('d', [20439.3, 16139.0, 13494.1, 10286.9, 7665.9, 5926.6, 4486.9, 3323.6, 2492.4, 1855.9, 1372.6, 1009.7, 732.9, 531.0, 381.5, 273.6, 195.6, 137.9, 95.8, 64.7, 46.4, 31.1, 21.7, 14.9, 9.7, 7.4, 5.0, 3.5, 2.7, 2.0, 1.4, 1.1, 0.8, 0.7, 0.7, 0.6, 0.6, 0.4])

#off['trackerMET'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0])
#offrate['trackerMET'] =  array('d', [13268.1, 8221.6, 4763.9, 2641.6, 1442.3, 792.3, 439.4, 255.7, 152.8, 96.2, 63.0, 44.5, 30.9, 23.1, 16.9, 13.5, 11.0, 8.9])
#onl['trackerMET'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0])
#onlrate['trackerMET'] =  array('d', [162.4, 22.8, 7.7, 4.4, 3.0, 2.4, 0.9, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1])

off['trackerMET'] =  array('d', [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0])
offrate['trackerMET'] =  array('d', [26259.8, 25030.5, 23749.2, 22389.2, 21047.6, 19697.7, 18342.4, 17018.3, 15722.9, 14473.4, 13268.1, 12134.5, 11055.6, 10047.1, 9103.1, 8221.6, 7395.4, 6641.8, 5944.8, 5328.9, 4763.9, 4244.5, 3768.4, 3345.3, 2970.2, 2641.6, 2342.7, 2075.6, 1844.6, 1630.8, 1442.3, 1273.2, 1128.9, 1005.2, 890.2, 792.3, 705.8, 623.2, 553.1, 494.3, 439.4, 393.5, 353.4, 316.0, 285.2, 255.7, 232.2, 209.2, 188.7, 167.7, 152.8, 140.1, 126.6, 117.0, 104.9, 96.2, 88.4, 81.3, 73.2, 67.9, 63.0, 58.1, 54.4, 51.2, 47.4, 44.5, 41.3, 38.6, 35.7, 33.3, 30.9, 28.5, 26.8, 25.4, 23.9, 23.1, 21.9, 20.2, 19.0, 17.8, 16.9, 15.9, 14.9, 14.5, 14.2, 13.5, 13.1, 12.6, 12.0, 11.1, 11.0, 10.5, 10.3, 9.6, 9.3, 8.9, 8.7, 8.5, 8.3, 8.1])
onl['trackerMET'] =  array('d', [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0])
onlrate['trackerMET'] =  array('d', [28284.6, 23248.7, 16255.6, 10240.4, 5930.4, 3253.3, 1743.2, 927.3, 500.2, 281.4, 162.4, 100.1, 64.0, 44.7, 30.4, 22.8, 16.1, 13.2, 10.5, 8.7, 7.7, 7.2, 6.1, 5.4, 4.7, 4.4, 4.2, 3.8, 3.5, 3.1, 3.0, 2.8, 2.5, 2.5, 2.4, 2.4, 2.2, 2.0, 1.9, 1.3, 0.9, 0.6, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

onl['trackerMET_FBE'] =  array('d', [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0])
onlrate['trackerMET_FBE'] =  array('d', [30655.7, 27369.0, 20398.5, 13074.0, 7426.3, 3881.7, 1965.9, 998.2, 522.0, 293.7, 174.0, 111.2, 76.2, 55.4, 40.8, 31.7, 25.8, 21.0, 17.7, 15.6, 14.1, 12.2, 10.8, 9.4, 8.2, 6.6, 3.3, 2.2, 1.7, 1.2, 1.0, 0.8, 0.7, 0.7, 0.6, 0.6, 0.6, 0.6, 0.5, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.2, 0.2, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])


off['trackerJet'] =  array('d', [40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0])
offrate['trackerJet'] =  array('d', [31036.5, 31036.5, 30915.7, 30680.9, 22388.0, 8428.3, 3884.0, 2065.0, 1190.5, 731.1, 474.2, 320.4, 221.5, 156.0, 112.2, 84.1, 64.2, 50.1, 38.1, 30.6])
onl['trackerJet'] =  array('d', [40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0])
onlrate['trackerJet'] =  array('d', [501.4, 118.3, 37.5, 15.1, 7.8, 5.0, 3.9, 3.0, 1.6, 0.5, 0.5, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3])

off['caloJet'] =  array('d', [40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0])
offrate['caloJet'] =  array('d', [31038.0, 31038.0, 26155.6, 5045.0, 948.3, 340.7, 159.3, 83.5, 47.2, 29.5, 18.6, 12.2, 7.5, 5.8, 3.8, 2.8, 2.1, 1.6, 1.2, 1.1])
onl['caloJet'] =  array('d', [40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0])
onlrate['caloJet'] =  array('d', [12842.3, 1955.9, 539.3, 209.7, 92.6, 46.7, 27.1, 15.5, 8.6, 5.4, 3.7, 2.3, 1.5, 1.2, 1.1, 0.8, 0.5, 0.5, 0.3, 0.1])

off['caloJetExt'] =  array('d', [40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0])
offrate['caloJetExt'] =  array('d', [31038.0, 31038.0, 26303.9, 5485.8, 1167.2, 449.1, 214.9, 112.4, 65.1, 39.4, 24.8, 16.1, 10.5, 7.7, 5.0, 3.5, 2.5, 2.1, 1.5, 1.4])
onl['caloJetExt'] =  array('d', [40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0])
onlrate['caloJetExt'] =  array('d', [14705.0, 2358.7, 643.6, 242.6, 105.7, 51.5, 29.5, 16.9, 8.9, 5.5, 4.0, 2.6, 1.5, 1.2, 1.1, 0.8, 0.5, 0.5, 0.3, 0.1])

off['NNPuppiTauLoose'] =  array('d', [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0])
offrate['NNPuppiTauLoose'] =  array('d', [2406.3, 2182.7, 1846.7, 1449.3, 1049.6, 726.3, 525.7, 387.3, 294.1, 224.5, 179.0, 141.8, 115.2, 93.5, 79.0, 69.1, 58.7, 50.9, 44.7, 39.4, 35.3, 31.0, 28.4, 25.8, 23.4, 21.3, 19.1, 17.6, 16.3, 15.4])
onl['NNPuppiTauLoose'] =  array('d', [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0])
onlrate['NNPuppiTauLoose'] =  array('d', [2402.4, 1934.4, 1278.2, 725.6, 438.7, 284.4, 190.6, 132.9, 96.0, 74.8, 57.1, 46.3, 38.3, 31.8, 27.7, 24.1, 19.9, 17.2, 14.9, 13.5, 12.0, 10.6, 9.5, 9.0, 8.3, 7.4, 6.5, 6.0, 5.3, 4.8])

off['CaloTau'] =  array('d', [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0])
offrate['CaloTau'] =  array('d', [31037.9, 31037.6, 30857.6, 27058.4, 17349.8, 9348.6, 5053.3, 2881.1, 1778.3, 1179.3, 826.7, 601.5, 450.0, 345.6, 271.4, 211.1, 171.7, 139.3, 114.4, 93.6, 78.7, 66.1, 56.7, 48.0, 41.4, 36.2, 31.0, 26.8, 22.7, 19.1])
onl['CaloTau'] =  array('d', [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0])
onlrate['CaloTau'] =  array('d', [31037.9, 30775.8, 23729.0, 11431.6, 5007.3, 2421.8, 1342.8, 822.4, 540.0, 375.5, 269.2, 195.9, 147.9, 113.5, 88.8, 70.5, 56.4, 46.5, 37.8, 31.4, 25.1, 20.5, 16.6, 14.6, 12.8, 11.1, 9.3, 7.9, 6.9, 6.4])

##these are inclusive ttbar
off['puppiPhase1HT'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0])
offrate['puppiPhase1HT'] =  array('d', [31038.0, 5536.7, 2769.4, 1507.1, 887.9, 518.6, 313.3, 193.6, 123.4, 84.0, 60.1, 42.3, 29.8, 21.5, 16.2, 12.7, 10.0, 8.1, 6.2, 4.7, 4.0, 3.2, 2.5, 2.0, 1.8, 1.4, 1.2, 1.1, 0.8, 0.6, 0.6, 0.6, 0.6, 0.5, 0.4, 0.4, 0.3, 0.3])
onl['puppiPhase1HT'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0])
onlrate['puppiPhase1HT'] =  array('d', [2133.6, 1254.3, 692.4, 395.6, 232.4, 143.8, 91.7, 63.8, 44.1, 30.2, 21.2, 15.8, 12.2, 9.6, 7.6, 5.7, 4.3, 3.5, 2.6, 2.0, 1.8, 1.4, 1.2, 1.1, 0.8, 0.6, 0.6, 0.6, 0.6, 0.5, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3])


off['trackerHT'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0])
offrate['trackerHT'] =  array('d', [29767.4, 21615.9, 12100.4, 6120.9, 3106.2, 1659.0, 940.2, 566.0, 358.2, 234.9, 164.3, 114.5, 83.6, 62.0, 47.0, 36.0, 28.2, 22.3, 18.8, 15.4, 12.9, 11.0, 8.9, 7.6, 6.6, 5.5, 4.6, 3.9, 2.9, 2.6, 2.2, 1.5, 1.4, 1.3, 1.2, 1.1, 1.1, 0.9])
onl['trackerHT'] =  array('d', [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0])
onlrate['trackerHT'] =  array('d', [2809.5, 586.8, 182.1, 73.1, 34.2, 18.8, 11.3, 7.0, 4.4, 2.6, 1.4, 1.1, 0.9, 0.6, 0.6, 0.5, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1])

