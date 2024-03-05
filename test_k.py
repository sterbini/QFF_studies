# %% Preparing the knobs for the matching
from cpymad.madx import Madx


my_mad = Madx()
my_mad.call('/afs/cern.ch/eng/lhc/optics/runIII/LHC_LS2_2021-07-02.seq')


my_list = [
'MQY.4R1.B2',
'MQML.5R1.B2',
'MQML.6R1.B2',
'MQM.A7R1.B2',
'MQM.B7R1.B2',
'MQML.8R1.B2',
'MQM.9R1.B2',
'MQML.10R1.B2',
'MQTLI.11R1.B2',
'MQT.12R1.B2',
'MQT.13R1.B2',
'MQY.4L1.B2',
'MQML.5L1.B2',
'MQML.6L1.B2',
'MQM.A7L1.B2',
'MQM.B7L1.B2',
'MQML.8L1.B2',
'MQM.9L1.B2',
'MQML.10L1.B2',
'MQTLI.11L1.B2',
'MQT.12L1.B2',
'MQT.13L1.B2',
]
# %%
for ii in my_list:
    print(f"{ii}:\t [{my_mad.elements[ii].kmin}, {my_mad.elements[ii].kmax}]")
# %%
