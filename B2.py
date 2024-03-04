# %%
#/afs/cern.ch/eng/lhc/optics/runIII/RunIII_dev/Proton_2024/README
# opticsfile.43
import json
from matplotlib import pyplot as plt
import xtrack as xt
import numpy as np

if False:
        from build_distr_and_collider import *
        build_distr_and_collider()

# %%
my_beam = 'b2' # this is the beam 4 in madx sense
my_optics = 30 # beta* at which the knob is matched

if my_optics==30:
        collider = xt.Multiline.from_json('./collider/collider_43.json')
elif my_optics==28:
        collider = xt.Multiline.from_json('./collider/collider_44.json')
elif my_optics==26:
        collider = xt.Multiline.from_json('./collider/collider_45.json')
elif my_optics==24:
        collider = xt.Multiline.from_json('./collider/collider_46.json')
elif my_optics==22:
        collider = xt.Multiline.from_json('./collider/collider_47.json')
elif my_optics==20:
        collider = xt.Multiline.from_json('./collider/collider_48.json')
elif my_optics==60:
        collider = xt.Multiline.from_json('./collider/collider_34.json')
else:
        raise ValueError('The optics is not available')

# remove the my_list the elements mbas2.1l1 and mbas2.1r1 from collider 
s_ip = collider[f'lhc{my_beam}'].get_s_position(at_elements='ip1')
line = collider[f'lhc{my_beam}']#.filter_elements(exclude_types_starting_with='Solenoid')
sol= line['mbas2.1l1'] 
line.element_dict['mbas2.1l1'] =xt.Drift(length=sol.length)
sol= line['mbas2.1r1'] 
line.element_dict['mbas2.1r1'] =xt.Drift(length=sol.length)

# %%
# Check of the tctpv center
# position of the wire in the tctpv.4r1.b2 (center of the TCTPV)
print('The wire in the TCTPV.4R1.B2 is at', line.get_s_position(f'tctpv.4r1.b2')-line.get_s_position(f'ip1'), 'm from the IP1.')
print('The wire in the TCTPH.4R5.B2 is at', line.get_s_position(f'tctph.4r5.b2')-line.get_s_position(f'ip5'), 'm from the IP5.')



# %%
# installation of the wires in the sequence

line.vars[f'i_wire_ip1.{my_beam}'] = 0.0 
line.vars[f'd_wire_ip1.{my_beam}'] = 0.01 

line.vars[f'i_wire_ip5.{my_beam}'] = 0.0
line.vars[f'd_wire_ip5.{my_beam}'] = 0.01

# upper wire in IR1
line.insert_element(name=f'bbwc.t.4r1.{my_beam}',
                    element=xt.Wire(
                        L_phy=1, 
                        L_int=2,
                        current= 0.0,
                        xma=0.0, 
                        yma= 1 # very far from the beam
                        ),
                    at=f'tctpv.4r1.{my_beam}')
# bottom wire in IR1
line.insert_element(name=f'bbwc.b.4r1.{my_beam}',
                    element=xt.Wire(
                        L_phy=1, 
                        L_int=2, 
                        current=0.0, 
                        xma=0.0, 
                        yma=-1  # very far from the beam
                        ),
                    at=f'tctpv.4r1.{my_beam}')

# upper wire in IR5
line.insert_element(name=f'bbwc.e.4r5.{my_beam}',
                    element=xt.Wire(
                        L_phy=1, 
                        L_int=2,
                        current= 0.0,
                        xma=1, # very far from the beam
                        yma=0.0 
                        ),
                    at=f'tctph.4r5.{my_beam}')
# bottom wire in IR5
line.insert_element(name=f'bbwc.i.4r5.{my_beam}',
                    element=xt.Wire(
                        L_phy=1, 
                        L_int=2, 
                        current=0.0, 
                        xma=-1, # very far from the beam
                        yma=0.0  
                        ),
                    at=f'tctph.4r5.{my_beam}')

s_ip = line.get_s_position(at_elements='ip1')

# I am adding markers to check the beta-beating close to the IP 
for ii in range(1,60):
        print(f'Adding marker at {ii/20} m from IP1 (right)')
        line.insert_element(name=f'marker_finer.{ii}.r1.{my_beam}',
                    element=xt.Marker(),
                    at_s=s_ip + ii/20)
for ii in range(1,60):
        print(f'Adding marker at {ii/20} m from IP1 (left)')
        line.insert_element(name=f'marker_finer.{ii}.l1.{my_beam}',
                    element=xt.Marker(),
                    at_s=s_ip - ii/20)


# %%
# 't' is for top, 'b' is for bottom, 'e' is for external, 'i' is for internal
# this are the wires for the beam 2        
line.vars['co_y_wire_ip1.b2'] = 0
line.vars['co_x_wire_ip1.b2'] = 0
line.vars['co_y_wire_ip5.b2'] = 0
line.vars['co_x_wire_ip5.b2'] = 0        

line.element_refs['bbwc.t.4r1.b2'].current = line.vars['i_wire_ip1.b2']
line.element_refs['bbwc.t.4r1.b2'].yma = line.vars['d_wire_ip1.b2'] + line.vars['co_y_wire_ip1.b2']

line.element_refs['bbwc.b.4r1.b2'].current = line.vars['i_wire_ip1.b2'] 
line.element_refs['bbwc.b.4r1.b2'].yma = -line.vars['d_wire_ip1.b2'] + line.vars['co_y_wire_ip1.b2']

line.element_refs['bbwc.e.4r5.b2'].current = line.vars['i_wire_ip5.b2']
line.element_refs['bbwc.e.4r5.b2'].xma = line.vars['d_wire_ip5.b2'] + line.vars['co_x_wire_ip5.b2']

line.element_refs['bbwc.i.4r5.b2'].current = line.vars['i_wire_ip5.b2']
line.element_refs['bbwc.i.4r5.b2'].xma = -line.vars['d_wire_ip5.b2'] + line.vars['co_x_wire_ip5.b2']


# %%
# flat orbit
for ii in [2,8]:
        collider.vars[f'on_x{ii}h'] = 0.0
        collider.vars[f'on_x{ii}v'] = 0.0
        collider.vars[f'on_sep{ii}h'] = 0.0
        collider.vars[f'on_sep{ii}v'] = 0.0
# %%
for ii in [1,2,5,8]:
        print(8*'*', f'IP{ii}', 8*'*')
        if ii in [2,8]:
                print(f'on_sep{ii}h:\t ', collider.vars[f'on_sep{ii}h']._get_value())
                print(f'on_sep{ii}v:\t ', collider.vars[f'on_sep{ii}v']._get_value())
        else:
                print(f'on_x{ii}:\t\t ', collider.vars[f'on_x{ii}']._get_value())
                print(f'on_sep{ii}:\t ', collider.vars[f'on_sep{ii}']._get_value())
        print(f'on_oh{ii}:\t\t ', collider.vars[f'on_oh{ii}']._get_value())
        print(f'on_ov{ii}:\t\t ', collider.vars[f'on_ov{ii}']._get_value())
        print(f'on_a{ii}:\t\t ', collider.vars[f'on_a{ii}']._get_value())

print(8*'*', 'others settings', 8*'*')
print('on_alice_normalized:\t', collider.vars['on_alice_normalized']._get_value())
print('on_lhcb_normalized:\t', collider.vars['on_lhcb_normalized']._get_value())
print('on_disp:\t\t', collider.vars['on_disp']._get_value())
# %%
tw_b2 = line.twiss(method='4d')
from matplotlib import pyplot as plt
plt.plot(tw_b2['s'], tw_b2['x'], label='x')
plt.plot(tw_b2['s'], tw_b2['y'], label='y')
# set the xticks only at the IPs
plt.xticks([tw_b2['s','ip1'], tw_b2['s','ip2'], tw_b2['s','ip5'], tw_b2['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])
plt.title('Beam 1')
plt.legend()
assert tw_b2['x'].std() < 1e-8
assert tw_b2['y'].std() < 1e-8
# %%
# set on_x1 to 160
collider.vars['on_x1'] = 160.0
tw_b2 = line.twiss(method='4d')
plt.plot(tw_b2['s']- tw_b2.rows['ip1']['s'], tw_b2['x'], label='x')
plt.plot(tw_b2['s']- tw_b2.rows['ip1']['s'], tw_b2['y'], label='y')

plt.axvline(x=tw_b2['s','bbwc.b.4r1.b2']- tw_b2.rows['ip1']['s'],
             color='r', 
             linestyle='-',
             label='bbwc.b.4r1.b2')
            
plt.axvline(x=tw_b2['s','mqml.5r1.b2']- tw_b2.rows['ip1']['s'],
             color='b', 
             linestyle='--', label='mqml.5r1.b2')

plt.axvline(x=tw_b2['s','mqml.5l1.b2']- tw_b2.rows['ip1']['s'],
             color='b', 
             linestyle='-.', label='mqml.5l1.b2')
plt.xlim(-300, 300)
#set back to 0
collider.vars['on_x1'] = 0.0
plt.xlabel('s [m]')
plt.ylabel('[m]')
plt.legend()

# %%
# set on_x1 to 160
collider.vars['on_x5'] = 160.0
tw_b2 = line.twiss(method='4d')
plt.plot(tw_b2['s']- tw_b2.rows['ip5']['s'], tw_b2['x'], label='x')
plt.plot(tw_b2['s']- tw_b2.rows['ip5']['s'], tw_b2['y'], label='y')

plt.axvline(x=tw_b2['s','bbwc.e.4r5.b2']- tw_b2.rows['ip5']['s'],
             color='r', 
             linestyle='-',
             label='bbwc.e.4r5.b2')
            

plt.xlim(-300, 300)
#set back to 0
collider.vars['on_x5'] = 0.0
plt.xlabel('s [m]')
plt.ylabel('[m]')
plt.legend()

# %%
collider.vars['on_sep1'] = 1.0
tw_b2 = line.twiss(method='4d')
plt.plot(tw_b2['s']- tw_b2.rows['ip1']['s'], tw_b2['x'], label='x')
plt.plot(tw_b2['s']- tw_b2.rows['ip1']['s'], tw_b2['y'], label='y')

plt.axvline(x=tw_b2['s','bbwc.b.4r1.b2']- tw_b2.rows['ip1']['s'],
             color='r', 
             linestyle='-',
             label='bbwc.b.4r1.b2')
            
plt.xlim(-300, 300)
collider.vars['on_sep1'] = 0.0
plt.xlabel('s [m]')
plt.ylabel('[m]')
plt.legend()

# %%
collider.vars['on_sep5'] = 1.0
tw_b2 = line.twiss(method='4d')
plt.plot(tw_b2['s']- tw_b2.rows['ip5']['s'], tw_b2['x'], label='x')
plt.plot(tw_b2['s']- tw_b2.rows['ip5']['s'], tw_b2['y'], label='y')

plt.axvline(x=tw_b2['s','bbwc.e.4r5.b2']- tw_b2.rows['ip5']['s'],
             color='r', 
             linestyle='-',
             label='bbwc.e.4r5.b2')
            
plt.xlim(-300, 300)
collider.vars['on_sep5'] = 0.0
plt.xlabel('s [m]')
plt.ylabel('[m]')
plt.legend()
# %%
for ip in [1,5]:
        assert np.isclose(tw_b2['betx',f'ip{ip}'], collider.vars[f'betx_ip{ip}']._get_value(), rtol=1e-4)
        assert np.isclose(tw_b2['bety',f'ip{ip}'], collider.vars[f'bety_ip{ip}']._get_value(), rtol=1e-4)
print('Beam 1')
assert np.isclose(tw_b2.qx, 62.31, atol=1e-6)
assert np.isclose(tw_b2.qy, 60.32, atol=1e-6)
tw_b2[['betx','alfx','bety','alfy','x','px','y','py'],'ip.*'].to_pandas()

# %%
# For reference
# these are the bbwc installed in the sequence
print('These are the bbwc installed in the sequence', [ii for ii in collider.lhcb2.element_names if ('bbwc' in ii) ])



# %%
# here we need to load the value of the collider at 30 cm to compute the beam-wire distance at 30 cm
collider_at_30cm =  xt.Multiline.from_json('./collider/collider_43.json')
for ii in [2,8]:
        collider_at_30cm.vars[f'on_x{ii}h'] = 0.0
        collider_at_30cm.vars[f'on_x{ii}v'] = 0.0
        collider_at_30cm.vars[f'on_sep{ii}h'] = 0.0
        collider_at_30cm.vars[f'on_sep{ii}v'] = 0.0

epsilon_collimation = 3.5e-6
wire_retraction = 0.003
tct_opening_in_sigma = 8 
proton_mass_in_GeV = 0.93827208816
# relativistic beta of a proton at 6.8 TeV
beta_rel_proton = np.sqrt(1- (proton_mass_in_GeV/collider.vars['nrj']._get_value())**2)
# compute gamma_rel_proton
gamma_rel_proton = 1/np.sqrt(1-beta_rel_proton**2)
tw_special = collider_at_30cm['lhcb2'].twiss(method='4d')

sigma_y_at_tctpv_4r1_b2 = np.sqrt(tw_special['bety','tctpv.4r1.b2']
                                * epsilon_collimation
                                / beta_rel_proton
                                / gamma_rel_proton)

sigma_x_at_tctph_4r5_b2 = np.sqrt(tw_special['betx','tctph.4r5.b2']
                                * epsilon_collimation
                                / beta_rel_proton
                                / gamma_rel_proton)

line.vars['d_wire_ip1.b2'] = tct_opening_in_sigma * sigma_y_at_tctpv_4r1_b2 + wire_retraction
line.vars['d_wire_ip5.b2'] = tct_opening_in_sigma * sigma_x_at_tctph_4r5_b2 + wire_retraction


print('The beam energy is', collider.vars['nrj']._get_value(), 'GeV.')
print('The TCT opening assumed is', tct_opening_in_sigma,'collimation sigma.')
print('The wire retraction is', wire_retraction*1e3, 'mm.')
print('The wire-to-beam distance in IP1 is', line.vars['d_wire_ip1.b2']._get_value()*1e3, 'mm.')
print('The wire-to-beam distance in IP5 is', line.vars['d_wire_ip5.b2']._get_value()*1e3, 'mm.')
# %%
# to test that there is not orbit effect
line.vars['i_wire_ip1.b2'] = 350.0
line.vars['i_wire_ip5.b2'] = 350.0

tw_b2 = line.twiss(method='4d')
# to check that the orbit is flat
assert tw_b2['x'].std() < 1e-14
assert tw_b2['y'].std() < 1e-14

print(tw_b2.qx, tw_b2.qy)

# %%
line.vars['i_wire_ip1.b2'] = 0.0
line.vars['i_wire_ip5.b2'] = 0.0
tw_b2 = line.twiss(method='4d')

qx0 = tw_b2.qx
qy0 = tw_b2.qy
betx_0 = tw_b2.betx
bety_0 = tw_b2.bety


i_range = np.linspace(0, 350, 10)
delta_qx = []
delta_qy = []
for i_wire in i_range:
        line.vars['i_wire_ip1.b2'] = i_wire
        tw_b2 = line.twiss(method='4d')
        delta_qx.append(tw_b2.qx-qx0)
        delta_qy.append(tw_b2.qy-qy0)
        print(delta_qx[-1], delta_qy[-1])
        assert tw_b2['x'].std() < 1e-8 
        assert tw_b2['y'].std() < 1e-8

plt.figure()
plt.plot(tw_b2.s, (tw_b2.betx-betx_0)/betx_0, label='$\Delta\\beta$x/$\\beta$x$_0$')
plt.plot(tw_b2.s, (tw_b2.bety-bety_0)/bety_0, label='$\Delta\\beta$y/$\\beta$y$_0$')
plt.legend()
plt.xlabel('s [m]')
plt.ylabel('relative $\\beta$-beating')
plt.title('Wires at IR1 for B2 (350 A, TCT at 8$\sigma$)')
plt.ylim(-.16, 0.16)
plt.grid()
# set ticks only at the IPs
plt.xticks([tw_b2['s','ip1'], tw_b2['s','ip2'], tw_b2['s','ip5'], tw_b2['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])

plt.figure()
plt.plot(i_range, delta_qx, 's-', label='$\Delta$Qx')
plt.plot(i_range, delta_qy, 's-', label='$\Delta$Qy')
plt.grid()
plt.legend()
plt.xlabel('Current [A]')
plt.ylabel('Tune shift [-]')
plt.title('Wires in IR1 for B2 (350 A, TCT at 8$\sigma$)')
plt.ylim(-.01, 0.021)

# %%
line.vars['i_wire_ip1.b2'] = 0.0
line.vars['i_wire_ip5.b2'] = 0.0
tw_b2 = line.twiss(method='4d')

qx0 = tw_b2.qx
qy0 = tw_b2.qy
betx_0 = tw_b2.betx
bety_0 = tw_b2.bety


i_range = np.linspace(0, 350, 10)
delta_qx = []
delta_qy = []
for i_wire in i_range:
        line.vars['i_wire_ip5.b2'] = i_wire
        tw_b2 = line.twiss(method='4d')
        delta_qx.append(tw_b2.qx-qx0)
        delta_qy.append(tw_b2.qy-qy0)
        print(delta_qx[-1], delta_qy[-1])
        assert tw_b2['x'].std() < 1e-8 
        assert tw_b2['y'].std() < 1e-8

plt.figure()
plt.plot(tw_b2.s, (tw_b2.betx-betx_0)/betx_0, label='$\Delta\\beta$x/$\\beta$x$_0$')
plt.plot(tw_b2.s, (tw_b2.bety-bety_0)/bety_0, label='$\Delta\\beta$y/$\\beta$y$_0$')
plt.legend()
plt.xlabel('s [m]')
plt.ylabel('relative $\\beta$-beating')
plt.title('Wires at IR5 for B2 (350 A, TCT at 8$\sigma$)')
plt.ylim(-.16, 0.16)
plt.grid()
# set ticks only at the IPs
plt.xticks([tw_b2['s','ip1'], tw_b2['s','ip2'], tw_b2['s','ip5'], tw_b2['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])

plt.figure()
plt.plot(i_range, delta_qx, 's-', label='$\Delta$Qx')
plt.plot(i_range, delta_qy, 's-', label='$\Delta$Qy')
plt.grid()
plt.legend()
plt.xlabel('Current [A]')
plt.ylabel('Tune shift [-]')
plt.title('Wires at IR5 for B2 (350 A, TCT at 8$\sigma$)')
plt.ylim(-.01, 0.021)

# %% Preparing the knobs for the matching
my_k_list = [
           'kq5.l1b2',
           'kq5.r1b2', 
           'kq6.l1b2', 
           'kq6.r1b2',
           'kq7.l1b2',
           'kq7.r1b2',  
           'kq8.l1b2',
           'kq8.r1b2', 
           'kq9.l1b2', 
           'kq9.r1b2',
           'kq10.l1b2',
           'kq10.r1b2', 
           'kqtl11.r1b2', 
           'kqt12.r1b2', 
           'kqt13.r1b2',
           'kq4.l5b2',
           'kq4.r5b2',
           'kq5.l5b2',
           'kq5.r5b2', 
           'kq6.l5b2', 
           'kq6.r5b2',
           'kq7.l5b2',
           'kq7.r5b2',  
           'kq8.l5b2',
           'kq8.r5b2', 
           'kq9.l5b2', 
           'kq9.r5b2',
           'kq10.l5b2',
           'kq10.r5b2', 
           'kqtl11.r5b2', 
           'kqt12.r5b2', 
           'kqt13.r5b2',
           ]

# from /afs/cern.ch/eng/lhc/optics/runIII/LHC_LS2_2021-07-02.seq
limits_dict = { 'kq5.l1b2': ['kmin_mqml_4.5k','kmax_mqml_4.5k'], # MQML modified
                'kq5.r1b2': ['kmin_mqml_4.5k','kmax_mqml_4.5k'], # MQML
                'kq6.l1b2': ['kmin_mqml_4.5k','kmax_mqml_4.5k'], # MQML
                'kq6.r1b2': ['kmin_mqml_4.5k','kmax_mqml_4.5k'], # MQML
                'kq7.l1b2': ['kmin_mqm','kmax_mqm'], # MQM
                'kq7.r1b2': ['kmin_mqm','kmax_mqm'], # MQM
                'kq8.l1b2': ['kmin_mqml','kmax_mqml'], # MQML
                'kq8.r1b2': ['kmin_mqml','kmax_mqml'], # MQML
                'kq9.l1b2': ['kmin_mqmc','kmax_mqmc'], # MQMC
                'kq9.r1b2': ['kmin_mqmc','kmax_mqmc'], # MQMC
                'kq10.l1b2': ['kmin_mqml', 'kmax_mqml'], # MQML
                'kq10.r1b2': ['kmin_mqml', 'kmax_mqml'], # MQML
                'kqtl11.r1b2': ['kmin_mqtli','kmax_mqtli'], # MQTLI
                'kqt12.r1b2': ['kmin_mqt','kmax_mqt'], # MQT
                'kqt13.r1b2': ['kmin_mqt','kmax_mqt'], # MQT
                'kq4.l5b2': ['kmin_mqy_4.5k','kmax_mqy_4.5k'], # MQY
                'kq4.r5b2': ['kmin_mqy_4.5k','kmax_mqy_4.5k'], # MQY
                'kq5.l5b2': ['kmin_mqml_4.5k','kmax_mqml_4.5k'], # MQML
                'kq5.r5b2': ['kmin_mqml_4.5k','kmax_mqml_4.5k'], # MQML
                'kq6.l5b2': ['kmin_mqml_4.5k','kmax_mqml_4.5k'], # MQML
                'kq6.r5b2': ['kmin_mqml_4.5k','kmax_mqml_4.5k'], # MQML
                'kq7.l5b2': ['kmin_mqm','kmax_mqm'], # MQM
                'kq7.r5b2': ['kmin_mqm','kmax_mqm'], # MQM
                'kq8.l5b2': ['kmin_mqml','kmax_mqml'], # MQML
                'kq8.r5b2': ['kmin_mqml','kmax_mqml'], # MQML
                'kq9.l5b2': ['kmin_mqmc','kmax_mqmc'], # MQMC
                'kq9.r5b2': ['kmin_mqmc','kmax_mqmc'], # MQMC
                'kq10.l5b2': ['kmin_mqml', 'kmax_mqml'], # MQML
                'kq10.r5b2': ['kmin_mqml', 'kmax_mqml'], # MQML
                'kqtl11.r5b2': ['kmin_mqtli','kmax_mqtli'], # MQTLI
                'kqt12.r5b2': ['kmin_mqt','kmax_mqt'], #MQT
                'kqt13.r5b2': ['kmin_mqt','kmax_mqt'], #MQT
               }

def print_k_summary(my_k_list):
        for ii in my_k_list:
                if not collider.vars[f'{ii}_delta']._get_value() == 0.0: 
                        print(ii, f'= {collider.vars[ii]._get_value()} =', collider.vars[f'{ii}_0']._get_value(),
                        ' + (', collider.vars[f'{ii}_delta']._get_value(), '), delta =',
                        collider.vars[f'{ii}_delta']._get_value()/collider.vars[f'{ii}_0']._get_value()*100, '%')
                        if (collider.vars[limits_dict[ii][0]]._get_value()/(3.3356*6800)<np.abs(collider.vars[ii]._get_value()) <
                        collider.vars[limits_dict[ii][1]]._get_value()/(3.3356*6800)):
                                print(f'{ii} is  within the limits',
                                       f'assuming [{collider.vars[limits_dict[ii][0]]._get_value()/(3.3356*6800)},',
                                       f'{collider.vars[limits_dict[ii][1]]._get_value()/(3.3356*6800)}]\n')
                        else:
                                print(ii, 'is OUTSIDE the limits\n')

def reset_delta_k(my_k_list):
        for kk in my_k_list:
                collider.vars[f'{kk}_delta'] = 0.000000

reset_delta_k(my_k_list)

for ii in my_k_list:
        collider.vars[f'{ii}_0'] = collider.vars[ii]._get_value()
        collider.vars[f'{ii}_delta'] = 0.000000
        if 'r1' in ii:
                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip1.b2']/350
        if 'l1' in ii:
                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip1.b2']/350
 
        if 'r5' in ii:
                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip5.b2']/350
        if 'l5' in ii:
                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip5.b2']/350


# %%

reset_delta_k(my_k_list)
line.vars['i_wire_ip1.b2'] = 0.0
line.vars['i_wire_ip5.b2'] = 0.0
tw_ref = line.twiss(method='4d')
def get_limits(my_k, factor=0.01):
        assert factor>0
        my_name = my_k.split('_delta')[0]
        nominal_value = line.vars[f'{my_name}_0']._get_value()
        if nominal_value>0:
                return (nominal_value*(-factor), nominal_value*(+factor))
        else:
                return (nominal_value*(+factor), nominal_value*(-factor))

for  match_ip1 in [True, False]:

        if match_ip1:
                line.vars['i_wire_ip1.b2'] = 350.0
                line.vars['i_wire_ip5.b2'] = 0.0

                tw_b2 = line.twiss(method='4d')
                print('Before matching')
                print(tw_b2.qx, tw_b2.qy)

                variables_dict = {
                                'kq5.r1b2_delta': {'limits': get_limits('kq5.r1b2_delta', 0.1),'step': 1e-8},
                                'kq6.r1b2_delta': {'limits': get_limits('kq6.r1b2_delta', 0.1),'step': 1e-8},
                                'kq7.r1b2_delta': {'limits': get_limits('kq7.r1b2_delta', 0.1),'step': 1e-8},
                                'kq8.r1b2_delta': {'limits': get_limits('kq8.r1b2_delta', 0.1),'step': 1e-8},
                                'kq9.r1b2_delta': {'limits': get_limits('kq9.r1b2_delta', 0.1),'step': 1e-8},
                                'kq10.r1b2_delta': {'limits': get_limits('kq10.r1b2_delta', 0.1),'step': 1e-8},
                                'kqtl11.r1b2_delta': {'limits': get_limits('kqtl11.r1b2_delta', 0.2),'step': 1e-8},
                                'kqt12.r1b2_delta': {'limits': get_limits('kqt12.r1b2_delta', 0.2),'step': 1e-8},
                                'kqt13.r1b2_delta': {'limits': get_limits('kqt13.r1b2_delta', 0.2),'step': 1e-8},
                                }       


                variables_list = [
                                'kq5.r1b2_delta', 
                                'kq6.r1b2_delta',
                                'kq7.r1b2_delta',
                                'kq8.r1b2_delta',
                                'kq9.r1b2_delta',
                                'kq10.r1b2_delta',
                                'kqtl11.r1b2_delta',
                                'kqt12.r1b2_delta',
                                'kqt13.r1b2_delta',                  
                                ]
                opt = line.match(
                solve=False, # <- prepare the match without running it
                start='e.ds.r1.b2', end='ip1',
                init=tw_ref, init_at=xt.START,
                method='4d',
                vary=
                        [xt.VaryList([ii], limits=variables_dict[ii]['limits'],  step=variables_dict[ii]['step']) for ii in variables_list ],
                targets=[
                        xt.TargetSet(['mux', 'muy','betx','bety','alfx','alfy',
                                'dx',
                                #'dpx',
                                ], value=tw_ref, at=xt.END),
                ])       
                my_ip = 1
                my_tct = 'tctpv.4r1.b2'
                my_i = 'i_wire_ip1.b2'

        else:
                line.vars['i_wire_ip1.b2'] = 0.0
                line.vars['i_wire_ip5.b2'] = 350.0
                tw_b2 = line.twiss(method='4d')
                print('Before matching')
                print(tw_b2.qx, tw_b2.qy)
        
                variables_dict = {
                                'kq4.r5b2_delta': {'limits': get_limits('kq4.r5b2_delta', 0.05),'step': 1e-8},
                                'kq5.r5b2_delta': {'limits': get_limits('kq5.r5b2_delta', 0.05),'step': 1e-8},
                                'kq6.r5b2_delta': {'limits': get_limits('kq6.r5b2_delta', 0.05),'step': 1e-8},
                                'kq7.r5b2_delta': {'limits': get_limits('kq7.r5b2_delta', 0.05),'step': 1e-8},
                                'kq8.r5b2_delta': {'limits': get_limits('kq8.r5b2_delta', 0.05),'step': 1e-8},
                                'kq9.r5b2_delta': {'limits': get_limits('kq9.r5b2_delta', 0.05),'step': 1e-8},
                                'kq10.r5b2_delta': {'limits': get_limits('kq10.r5b2_delta', 0.05),'step': 1e-8},
                                'kqtl11.r5b2_delta': {'limits': get_limits('kqtl11.r5b2_delta', 0.1),'step': 1e-8},
                                'kqt12.r5b2_delta': {'limits': get_limits('kqt12.r5b2_delta', 0.1),'step': 1e-8},
                                'kqt13.r5b2_delta': {'limits': get_limits('kqt13.r5b2_delta', 0.1),'step': 1e-8},
                                }      


                variables_list = ['kq4.r5b2_delta', 
                                'kq5.r5b2_delta', 
                                #'kq6.r5b2_delta',
                                'kq7.r5b2_delta',
                                'kq8.r5b2_delta',
                                'kq9.r5b2_delta',
                                'kq10.r5b2_delta',
                                'kqtl11.r5b2_delta',
                                'kqt12.r5b2_delta',
                                'kqt13.r5b2_delta',                  
                                ]
                
                opt = line.match(
                solve=False, # <- prepare the match without running it
                start='e.ds.r5.b2', end='ip5',
                init=tw_ref, init_at=xt.START,
                method='4d',
                vary=
                        [xt.VaryList([ii], limits=variables_dict[ii]['limits'],  step=variables_dict[ii]['step']) for ii in variables_list ],
                targets=[
                        xt.TargetSet(['mux', 'muy',
                                'betx','bety',
                                'alfx','alfy',
                                'dx',
                                'dpx',
                                ], value=tw_ref, at=xt.END),
                ])
                my_ip = 5
                my_tct = 'tctph.4r5.b2'
                my_i = 'i_wire_ip5.b2'
 
        opt.step(10)

        delta_dict = {}
        k0_dict = {}
        k_dict = {}
        kmin = {}
        kmax = {}
        k_relative_variation_percent = {}
        percentage_of_kmax = {}

        for ii in variables_list:
                name = ii.split('_delta')[0]
                k0_name = name+'_0' 
                delta_dict[name] = line.vars[ii]._get_value()

                k0_dict[name] = line.vars[k0_name]._get_value()
                k_dict[name] = k0_dict[name]+delta_dict[name]
                low_limit = line.vars[limits_dict[name][0]]._get_value()/6800/3.3356
                high_limit = line.vars[limits_dict[name][1]]._get_value()/6800/3.3356
                kmin[name] = low_limit
                kmax[name] = high_limit
                assert  kmin[name]<np.abs(k_dict[name])< kmax[name]
                k_relative_variation_percent[name] = delta_dict[name]/k0_dict[name]*100
                percentage_of_kmax[name] = np.abs(k_dict[name])/kmax[name]*100
        knob_dict = {}
        knob_dict['k_delta'] = delta_dict
        
        knob_dict['my_ip'] = my_ip
        knob_dict['my_optics'] = my_optics
        knob_dict['my_beam'] = my_beam
        knob_dict['i_wire_ip1.b2'] = line.vars['i_wire_ip1.b2']._get_value()
        knob_dict['i_wire_ip5.b2'] = line.vars['i_wire_ip5.b2']._get_value()
        knob_dict['tct_opening_in_sigma'] = tct_opening_in_sigma
        knob_dict['sigma_y_at_tctpv_4r1_b2'] = sigma_y_at_tctpv_4r1_b2
        knob_dict['sigma_x_at_tctph_4r5_b2'] = sigma_x_at_tctph_4r5_b2
        knob_dict['wire_retraction'] = wire_retraction
        knob_dict['d_wire_ip1.b2'] = line.vars['d_wire_ip1.b2']._get_value()
        knob_dict['d_wire_ip5.b2'] = line.vars['d_wire_ip5.b2']._get_value()

        knob_dict['k_0'] = k0_dict
        knob_dict['k'] = k_dict
        knob_dict['kmin'] = kmin
        knob_dict['kmax'] = kmax
        knob_dict['k_relative_variation_percent'] = k_relative_variation_percent
        knob_dict['percentage_of_kmax'] = percentage_of_kmax

        with open(f'knob_dict_350A_8sigma@30cm_ip{my_ip}_beta{my_optics}_{my_beam}.json', 'w') as f:
                json.dump(knob_dict, f, indent=4)


        tw_b2 = line.twiss(method='4d')
        print('After matching')
        print(tw_b2.qx, tw_b2.qy)
        
        plt.figure()
        plt.plot(tw_b2.s, (tw_b2.betx-tw_ref.betx)/tw_ref.betx, label='$\Delta\\beta_x/\\beta_{x0}$')
        plt.plot(tw_b2.s, (tw_b2.bety-tw_ref.bety)/tw_ref.bety, label='$\Delta\\beta_y\\beta_{y0}$')

        # set xticks only at the IPs
        plt.xticks([tw_b2['s','ip1'], tw_b2['s','ip2'], tw_b2['s','ip5'], tw_b2['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])
        plt.xlabel('s along the ring [no units]')
        plt.ylabel('relative $\\beta$-beating')
        plt.title(f'Wire at IP{my_ip}@350 A, {my_tct}@8$\sigma@30cm$ with optics at {my_optics} cm')
        plt.grid()
        plt.legend(loc='lower right')

        for ii in [my_ip]:
                plt.figure()
                s0 = tw_b2['s',f'ip{ii}']
                plt.plot(tw_b2.s-s0, (tw_b2.betx-tw_ref.betx)/tw_ref.betx, '.-', label='$\Delta\\beta_x/\\beta_{x0}$')
                plt.plot(tw_b2.s-s0, (tw_b2.bety-tw_ref.bety)/tw_ref.bety, '.-', label='$\Delta\\beta_y\\beta_{y0}$')
                plt.xlim(-700,5)
                plt.grid()
                plt.xlabel(f's from IP{my_ip} [m]')
                plt.ylabel('relative $\\beta$-beating')
                plt.title(f'Wire at IP{my_ip}@350 A, {my_tct}@8$\sigma@30cm$ with optics at {my_optics} cm')
                plt.legend(loc='upper left')

                

                # set the elements_list taking the list line.element_names from 's.ds.l1.b1' to 'ip1'
                element_list = line.element_names[line.element_names.index(f'e.ds.r{my_ip}.b2'):line.element_names.index(f'ip{my_ip}')+1]

                s_ref = line.get_s_position(f'ip{my_ip}')

                for element in element_list:
                        # if line.element_dict[element] belongs to xt.Quadrupole
                        if line[element].__class__.__name__ == 'Quadrupole':
                                my_dict = line[element].to_dict()
                                my_s = line.get_s_position(element)-s_ref
                                my_length = my_dict['length']
                                my_k = my_dict['k1']
                                # plot a red filled rectangle starting at my_s, 0 with length my_length and height my_k 
                                plt.gca().add_patch(plt.Rectangle((my_s, 0), my_length, my_k, color='r', alpha=0.5))
                        

                # plot betx and bety on the second vertical axis
                plt.twinx()
                plt.plot(tw_b2.s-s0, tw_b2.betx, 'r', label='$\\beta_x$')
                plt.plot(tw_b2.s-s0, tw_b2.bety, 'b', label='$\\beta_y$')
                plt.ylabel('$\\beta_x$ and $\\beta_y$ [m]')
                plt.legend(loc='lower right')

                plt.figure()
                plt.plot(tw_b2.s-s0, tw_b2.dx, '.-b', label='dx')
                plt.plot(tw_b2.s-s0, tw_b2.dy, '.-r', label='dy')
                plt.plot(tw_b2.s-s0, tw_ref.dx, 'ob', label='dx$_0$')
                plt.plot(tw_b2.s-s0, tw_ref.dy, 'or', label='dy$_0$')
                plt.xlim(-700,5)
                plt.ylim(-2.5, 0.5)
                plt.legend()
                plt.grid()
                plt.xlabel(f's from IP{my_ip} [m]')
                plt.ylabel('[m]')
                plt.title(f'Wire at IP{my_ip}@350 A, {my_tct}@8$\sigma@30cm$ with optics at {my_optics} cm')
                
                # print the k in my_k_list
                

        print(77*'*', f'\nFor {line.vars[my_i]._get_value()} A and 8 sigma at {my_tct}')
        print_k_summary(my_k_list)

# Check the linerarity of the tune knob


# %%
my_current = 50

if match_ip1:
        line.vars['i_wire_ip1.b2'] = my_current
        line.vars['i_wire_ip5.b2'] = 0.0
else:
        line.vars['i_wire_ip1.b2'] = 0.0
        line.vars['i_wire_ip5.b2'] = my_current

tw_b2 = line.twiss(method='4d')

print(tw_b2.qx-tw_ref.qx, tw_b2.qy-tw_ref.qy)
for ii in [my_ip]:
        plt.figure()
        s0 = tw_b2['s',f'ip{ii}']
        plt.plot(tw_b2.s-s0, (tw_b2.betx-tw_ref.betx)/tw_ref.betx, '.-', label='$\Delta\\beta_x/\\beta_{x0}$')
        plt.plot(tw_b2.s-s0, (tw_b2.bety-tw_ref.bety)/tw_ref.bety, '.-', label='$\Delta\\beta_y\\beta_{y0}$')
        #plt.xlim(-700,5)
        plt.grid()
        plt.xlabel(f's from IP{my_ip} [m]')
        plt.ylabel('relative $\\beta$-beating')
        plt.title(f'Wire at IP{my_ip}@{my_current} A, {my_tct}@8$\sigma@30cm$ with optics at {my_optics} cm')
        plt.legend(loc='upper left')
        plt.ylim(-.028, 0.028)

# %%
epsilon_geometric = epsilon_collimation/beta_rel_proton/gamma_rel_proton

line.vars['i_wire_ip1.b2'] = 0
line.vars['i_wire_ip5.b2'] = 0
line.vars['co_x_wire_ip1.b2'] = 0
line.vars['co_y_wire_ip1.b2'] = 0
line.vars['co_x_wire_ip5.b2'] = 0
line.vars['co_y_wire_ip5.b2'] = 0

line.vars['on_x1'] = 160.0

tw_b2 = line.twiss(method='4d')

line.vars['co_x_wire_ip1.b2'] = tw_b2['x', 'tctpv.4r1.b2']
line.vars['co_y_wire_ip1.b2'] = tw_b2['y', 'tctpv.4r1.b2']
line.vars['co_x_wire_ip5.b2'] = tw_b2['x', 'tctph.4r5.b2']
line.vars['co_y_wire_ip5.b2'] = tw_b2['y', 'tctph.4r5.b2']

line.vars['i_wire_ip1.b2'] = 350.0
tw_b2_wire_ip1_on = line.twiss(method='4d')
plt.plot(tw_b2_wire_ip1_on['mux'], (tw_b2_wire_ip1_on['x']- tw_b2['x'])/np.sqrt(tw_b2['betx']*epsilon_geometric), label='x')

plt.plot(tw_b2_wire_ip1_on['muy'], (tw_b2_wire_ip1_on['y']- tw_b2['y'])/np.sqrt(tw_b2['bety']*epsilon_geometric), label='y')
plt.ylabel('[$\sigma_{coll}$]')
plt.legend()
# show ticks of the IPs only
#plt.xticks([tw_b2_wire_ip1_on['s','ip1'], tw_b2_wire_ip1_on['s','ip5']], ['IP1', 'IP5'])
plt.title('Wires at IR1 for B2 (350 A, TCT at 8$\sigma$) and on_x1 = 160')
plt.xlabel('x and y phase [2$\pi$]')

line.vars['on_x1'] = 0.0
line.vars['on_x5'] = 0.0
line.vars['co_x_wire_ip1.b2'] = 0
line.vars['co_y_wire_ip1.b2'] = 0
line.vars['co_x_wire_ip5.b2'] = 0
line.vars['co_y_wire_ip5.b2'] = 0



# %%

line.vars['i_wire_ip1.b2'] = 0
line.vars['i_wire_ip5.b2'] = 0
line.vars['co_x_wire_ip1.b2'] = 0
line.vars['co_y_wire_ip1.b2'] = 0
line.vars['co_x_wire_ip5.b2'] = 0
line.vars['co_y_wire_ip5.b2'] = 0

line.vars['on_x5'] = 160.0
tw_b2 = line.twiss(method='4d')

line.vars['co_x_wire_ip1.b2'] = tw_b2['x', 'tctpv.4r1.b2']
line.vars['co_y_wire_ip1.b2'] = tw_b2['y', 'tctpv.4r1.b2']
line.vars['co_x_wire_ip5.b2'] = tw_b2['x', 'tctph.4r5.b2']
line.vars['co_y_wire_ip5.b2'] = tw_b2['y', 'tctph.4r5.b2']

line.vars['i_wire_ip5.b2'] = 350.0
tw_b2_wire_ip5_on = line.twiss(method='4d')

plt.plot(tw_b2_wire_ip5_on['mux'], (tw_b2_wire_ip5_on['x']- tw_b2['x'])/np.sqrt(tw_b2['betx']*epsilon_geometric), label='x')
plt.plot(tw_b2_wire_ip5_on['muy'], (tw_b2_wire_ip5_on['y']- tw_b2['y'])/np.sqrt(tw_b2['bety']*epsilon_geometric), label='y')
plt.ylabel('[$\sigma_{coll}$]')
plt.xlabel('x and y phase [2$\pi$]')
plt.legend()
# show ticks of the IPs only
#plt.xticks([tw_b2_wire_ip5_on['mux','ip1'], tw_b2_wire_ip5_on['s','ip5']], ['IP1', 'IP5'])
plt.title('Wires at IR1 for B2 (350 A, TCT at 8$\sigma$) and on_x5 = 160')

line.vars['on_x1'] = 0.0
line.vars['on_x5'] = 0.0
line.vars['co_x_wire_ip1.b2'] = 0
line.vars['co_y_wire_ip1.b2'] = 0
line.vars['co_x_wire_ip5.b2'] = 0
line.vars['co_y_wire_ip5.b2'] = 0
# %%
