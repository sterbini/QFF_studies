# %%
#/afs/cern.ch/eng/lhc/optics/runIII/RunIII_dev/Proton_2024/README
# opticsfile.43
from matplotlib import pyplot as plt
import xtrack as xt
import numpy as np
# import xobjects as xo

if False:
        from build_distr_and_collider import *
        build_distr_and_collider()

# %%
collider = xt.Multiline.from_json('./collider/collider_43.json')

my_beam = 'b2'
my_optics = 30

# remove the my_list the elements mbas2.1l1 and mbas2.1r1 from collider 
s_ip = collider[f'lhc{my_beam}'].get_s_position(at_elements='ip1')
line = collider[f'lhc{my_beam}']#.filter_elements(exclude_types_starting_with='Solenoid')
sol= line['mbas2.1l1'] 
line.element_dict['mbas2.1l1'] =xt.Drift(length=sol.length)
sol= line['mbas2.1r1'] 
line.element_dict['mbas2.1r1'] =xt.Drift(length=sol.length)
#line.build_tracker()


# %%
#line.discard_tracker()

line.vars[f'i_wire_ip1.{my_beam}'] = 0.0 
line.vars[f'd_wire_ip1.{my_beam}'] = 0.01 

line.vars[f'i_wire_ip5.{my_beam}'] = 0.0
line.vars[f'd_wire_ip5.{my_beam}'] = 0.01

# upper wire in IR1
line.insert_element(name=f'bbwc.t.4l1.{my_beam}',
                    element=xt.Wire(
                        L_phy=1, 
                        L_int=2,
                        current= 0.0,
                        xma=0.0, 
                        yma= 1 # very far from the beam
                        ),
                    at=f'tctpv.4l1.{my_beam}')
# bottom wire in IR1
line.insert_element(name=f'bbwc.b.4l1.{my_beam}',
                    element=xt.Wire(
                        L_phy=1, 
                        L_int=2, 
                        current=0.0, 
                        xma=0.0, 
                        yma=-1  # very far from the beam
                        ),
                    at=f'tctpv.4l1.{my_beam}')

# upper wire in IR5
line.insert_element(name=f'bbwc.e.4l5.{my_beam}',
                    element=xt.Wire(
                        L_phy=1, 
                        L_int=2,
                        current= 0.0,
                        xma=1, # very far from the beam
                        yma=0.0 
                        ),
                    at=f'tctph.4l5.{my_beam}')
# bottom wire in IR5
line.insert_element(name=f'bbwc.i.4l5.{my_beam}',
                    element=xt.Wire(
                        L_phy=1, 
                        L_int=2, 
                        current=0.0, 
                        xma=-1, # very far from the beam
                        yma=0.0  
                        ),
                    at=f'tctph.4l5.{my_beam}')

s_ip = line.get_s_position(at_elements='ip1')

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
#ctx = xo.ContextCpu(omp_num_threads='auto')
#line.build_tracker(_context=ctx)

# %%
# t is for top, b is for bottom, e is for external, i is for internal
# this are the wires for the beam 1        
line.element_refs['bbwc.t.4l1.b1'].current = line.vars['i_wire_ip1.b1']
line.element_refs['bbwc.t.4l1.b1'].yma = line.vars['d_wire_ip1.b1']

line.element_refs['bbwc.b.4l1.b1'].current = line.vars['i_wire_ip1.b1']
line.element_refs['bbwc.b.4l1.b1'].yma = -line.vars['d_wire_ip1.b1']

line.element_refs['bbwc.e.4l5.b1'].current = line.vars['i_wire_ip5.b1']
line.element_refs['bbwc.e.4l5.b1'].xma = line.vars['d_wire_ip5.b1']

line.element_refs['bbwc.i.4l5.b1'].current = line.vars['i_wire_ip5.b1']
line.element_refs['bbwc.i.4l5.b1'].xma = -line.vars['d_wire_ip5.b1']


# %%
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
tw_b1 = line.twiss(method='4d')

# %%
# for tw_b1
from matplotlib import pyplot as plt
plt.plot(tw_b1['s'], tw_b1['x'], label='x')
plt.plot(tw_b1['s'], tw_b1['y'], label='y')
# set the xticks only at the IPs
plt.xticks([tw_b1['s','ip1'], tw_b1['s','ip2'], tw_b1['s','ip5'], tw_b1['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])
plt.title('Beam 1')
plt.legend()
assert tw_b1['x'].std() < 1e-8
assert tw_b1['y'].std() < 1e-8
# %%
# set on_x1 to 160
collider.vars['on_x1'] = 160.0
tw_b1 = line.twiss(method='4d')
plt.plot(tw_b1['s']- tw_b1.rows['ip1']['s'], tw_b1['x'], label='x')
plt.plot(tw_b1['s']- tw_b1.rows['ip1']['s'], tw_b1['y'], label='y')

plt.axvline(x=tw_b1['s','bbwc.b.4l1.b1']- tw_b1.rows['ip1']['s'],
             color='r', 
             linestyle='-',
             label='bbwc.b.4l1.b1')
            
plt.axvline(x=tw_b1['s','mqml.5l1.b1']- tw_b1.rows['ip1']['s'],
             color='b', 
             linestyle='--')

plt.axvline(x=tw_b1['s','mqml.5r1.b1']- tw_b1.rows['ip1']['s'],
             color='b', 
             linestyle='--')
plt.xlim(-300, 300)
#set back to 0
collider.vars['on_x1'] = 0.0

# %%
collider.vars['on_sep1'] = 1.0
tw_b1 = line.twiss(method='4d')
plt.plot(tw_b1['s']- tw_b1.rows['ip1']['s'], tw_b1['x'], label='x')
plt.plot(tw_b1['s']- tw_b1.rows['ip1']['s'], tw_b1['y'], label='y')

plt.axvline(x=tw_b1['s','bbwc.b.4l1.b1']- tw_b1.rows['ip1']['s'],
             color='r', 
             linestyle='-',
             label='bbwc.b.4l1.b1')
            
plt.axvline(x=tw_b1['s','mqml.5l1.b1']- tw_b1.rows['ip1']['s'],
             color='b', 
             linestyle='--')

plt.axvline(x=tw_b1['s','mqml.5r1.b1']- tw_b1.rows['ip1']['s'],
             color='b', 
             linestyle='--')
plt.xlim(-300, 300)
#set back to 0
collider.vars['on_sep1'] = 0.0


# %%
for ip in [1,5]:
        assert np.isclose(tw_b1['betx',f'ip{ip}'], collider.vars[f'betx_ip{ip}']._get_value(), rtol=1e-4)
        assert np.isclose(tw_b1['bety',f'ip{ip}'], collider.vars[f'bety_ip{ip}']._get_value(), rtol=1e-4)
print('Beam 1')
assert np.isclose(tw_b1.qx, 62.31, atol=1e-6)
assert np.isclose(tw_b1.qy, 60.32, atol=1e-6)
# ask Gianni why the show is not working
tw_b1[['betx','alfx','bety','alfy','x','px','y','py'],'ip.*'].to_pandas()


# %%
# For reference
[ii for ii in collider.lhcb1.element_names if ('mqy.4r1.b1' in ii) ]


# %%
plt.plot(tw_b1.rows['e.ds.l1.b1':'s.ds.r1.b1']['s'] - tw_b1.rows['ip1']['s'],
        tw_b1.rows['e.ds.l1.b1':'s.ds.r1.b1']['betx'], label='betx')

plt.plot(tw_b1.rows['e.ds.l1.b1':'s.ds.r1.b1']['s']- tw_b1.rows['ip1']['s'],
        tw_b1.rows['e.ds.l1.b1':'s.ds.r1.b1']['bety'], label='bety')

# vertical line at the IP
plt.axvline(x=tw_b1['s','tctpv.4l1.b1']- tw_b1.rows['ip1']['s'],
             color='r', 
             linestyle='--')
# 
plt.axvline(x=tw_b1['s','bbwc.b.4l1.b1']- tw_b1.rows['ip1']['s'],
             color='r', 
             linestyle='-',
             label='tctpv.4l1.b1')

plt.axvline(x=tw_b1['s','mqy.4l1.b1']- tw_b1.rows['ip1']['s'],
             color='m', 
             linestyle='--',
             label='mqy.4l1.b1')

plt.axvline(x=tw_b1['s','mqy.4r1.b1']- tw_b1.rows['ip1']['s'],
             #dark green
             color='m', 
             linestyle='-.',
             label='mqy.4r1.b1')
            
plt.axvline(x=tw_b1['s','mqml.5l1.b1']- tw_b1.rows['ip1']['s'],
             color='b', 
             linestyle='--',
             label='mqml.5l1.b1')

plt.axvline(x=tw_b1['s','mqml.5r1.b1']- tw_b1.rows['ip1']['s'],
             color='b', 
             linestyle='-.',
             label='mqml.5r1.b1')

plt.axvline(x=tw_b1['s','mqml.6l1.b1']- tw_b1.rows['ip1']['s'],
             color='k', 
             linestyle='--',
             label='mqml.6l1.b1')

plt.axvline(x=tw_b1['s','mqml.6r1.b1']- tw_b1.rows['ip1']['s'],
             color='k', 
             linestyle='-.',
             label='mqml.6r1.b1')


plt.legend()
# %%
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
tw_special = collider_at_30cm['lhcb1'].twiss(method='4d')

sigma_y_at_tctpv_4l1_b1 = np.sqrt(tw_special['bety','tctpv.4l1.b1']
                                * epsilon_collimation
                                / beta_rel_proton
                                / gamma_rel_proton)

sigma_x_at_tctph_4l5_b1 = np.sqrt(tw_special['betx','tctph.4l5.b1']
                                * epsilon_collimation
                                / beta_rel_proton
                                / gamma_rel_proton)

line.vars['d_wire_ip1.b1'] = tct_opening_in_sigma * sigma_y_at_tctpv_4l1_b1 + wire_retraction
line.vars['d_wire_ip5.b1'] = tct_opening_in_sigma * sigma_x_at_tctph_4l5_b1 + wire_retraction


print('The beam energy is', collider.vars['nrj']._get_value(), 'GeV.')
print('The TCT opening assumed is', tct_opening_in_sigma,'collimation sigma.')
print('The wire retraction is', wire_retraction*1e3, 'mm.')
print('The wire-to-beam distance in IP1 is', line.vars['d_wire_ip1.b1']._get_value()*1e3, 'mm.')
print('The wire-to-beam distance in IP5 is', line.vars['d_wire_ip5.b1']._get_value()*1e3, 'mm.')
# %%
# to test that there is not orbit effect
line.vars['i_wire_ip1.b1'] = 350.0
line.vars['i_wire_ip5.b1'] = 350.0

tw_b1 = line.twiss(method='4d')
assert tw_b1['x'].std() < 1e-14
assert tw_b1['y'].std() < 1e-14

print(tw_b1.qx, tw_b1.qy)

# %%
line.vars['i_wire_ip1.b1'] = 0.0
line.vars['i_wire_ip5.b1'] = 0.0
tw_b1 = line.twiss(method='4d')

qx0 = tw_b1.qx
qy0 = tw_b1.qy
betx_0 = tw_b1.betx
bety_0 = tw_b1.bety


i_range = np.linspace(0, 350, 10)
delta_qx = []
delta_qy = []
for i_wire in i_range:
        line.vars['i_wire_ip1.b1'] = i_wire
        tw_b1 = line.twiss(method='4d')
        delta_qx.append(tw_b1.qx-qx0)
        delta_qy.append(tw_b1.qy-qy0)
        print(delta_qx[-1], delta_qy[-1])
        assert tw_b1['x'].std() < 1e-8 
        assert tw_b1['y'].std() < 1e-8

plt.figure()
plt.plot(tw_b1.s, (tw_b1.betx-betx_0)/betx_0, label='$\Delta\\beta$x/$\\beta$x$_0$')
plt.plot(tw_b1.s, (tw_b1.bety-bety_0)/bety_0, label='$\Delta\\beta$y/$\\beta$y$_0$')
plt.legend()
plt.xlabel('s [m]')
plt.ylabel('relative $\\beta$-beating')
plt.title('Effect of the wires at IP1 for B1 (350~A, 8$\sigma$)')
plt.ylim(-.16, 0.16)
plt.grid()

plt.figure()
plt.plot(i_range, delta_qx, label='$\Delta$Qx')
plt.plot(i_range, delta_qy, label='$\Delta$Qy')
plt.grid()
plt.legend()
plt.xlabel('Current [A]')
plt.ylabel('Tune shift [-]')
plt.title('Tune shift due to the wires at IP1 for B1')
plt.ylim(-.01, 0.021)

# %%
line.vars['i_wire_ip1.b1'] = 0.0
line.vars['i_wire_ip5.b1'] = 0.0
tw_b1 = line.twiss(method='4d')

qx0 = tw_b1.qx
qy0 = tw_b1.qy
betx_0 = tw_b1.betx
bety_0 = tw_b1.bety


i_range = np.linspace(0, 350, 10)
delta_qx = []
delta_qy = []
for i_wire in i_range:
        line.vars['i_wire_ip5.b1'] = i_wire
        tw_b1 = line.twiss(method='4d')
        delta_qx.append(tw_b1.qx-qx0)
        delta_qy.append(tw_b1.qy-qy0)
        print(delta_qx[-1], delta_qy[-1])
        assert tw_b1['x'].std() < 1e-8 
        assert tw_b1['y'].std() < 1e-8

plt.figure()
plt.plot(tw_b1.s, (tw_b1.betx-betx_0)/betx_0, label='$\Delta\\beta$x/$\\beta$x$_0$')
plt.plot(tw_b1.s, (tw_b1.bety-bety_0)/bety_0, label='$\Delta\\beta$y/$\\beta$y$_0$')
plt.legend()
plt.xlabel('s [m]')
plt.ylabel('relative $\\beta$-beating')
plt.title('Effect of the wires at IP5 for B1 (350 A, 8$\sigma$)')
plt.ylim(-.16, 0.16)
plt.grid()

plt.figure()
plt.plot(i_range, delta_qx, label='$\Delta$Qx')
plt.plot(i_range, delta_qy, label='$\Delta$Qy')
plt.grid()
plt.legend()
plt.xlabel('Current [A]')
plt.ylabel('Tune shift [-]')
plt.title('Tune shift due to the wires at IP5 for B1')
plt.ylim(-.01, 0.021)

# %%
my_k_list = [
           'kq5.l1b1',
           'kq5.r1b1', 
           'kq6.l1b1', 
           'kq6.r1b1',
           'kq7.l1b1',
           'kq7.r1b1',  
           'kq8.l1b1',
           'kq8.r1b1', 
           'kq9.l1b1', 
           'kq9.r1b1',
           'kq10.l1b1',
           'kq10.r1b1', 
           'kqtl11.l1b1', 
           'kqt12.l1b1', 
           'kqt13.l1b1',
           'kq4.l5b1',
           'kq4.r5b1',
           'kq5.l5b1',
           'kq5.r5b1', 
           'kq6.l5b1', 
           'kq6.r5b1',
           'kq7.l5b1',
           'kq7.r5b1',  
           'kq8.l5b1',
           'kq8.r5b1', 
           'kq9.l5b1', 
           'kq9.r5b1',
           'kq10.l5b1',
           'kq10.r5b1', 
           'kqtl11.l5b1', 
           'kqt12.l5b1', 
           'kqt13.l5b1',
           ]

limits_dict = { 'kq5.l1b1': ['kmin_mqml_4.5k','kmax_mqml_4.5k'],
                'kq5.r1b1': ['kmin_mqml_4.5k','kmax_mqml_4.5k'],
                'kq6.l1b1': ['kmin_mqml_4.5k','kmax_mqml_4.5k'],
                'kq6.r1b1': ['kmin_mqml_4.5k','kmax_mqml_4.5k'],
                'kq7.l1b1': ['kmin_mqm_4.5k','kmax_mqm_4.5k'],
                'kq7.r1b1': ['kmin_mqm_4.5k','kmax_mqm_4.5k'],
                'kq8.l1b1': ['kmin_mqml','kmax_mqml'],
                'kq8.r1b1': ['kmin_mqml','kmax_mqml'],
                'kq9.l1b1': ['kmin_mqml','kmax_mqml'],
                'kq9.r1b1': ['kmin_mqml','kmax_mqml'],
                'kq10.l1b1': ['kmin_mqml', 'kmax_mqml'], 
                'kq10.r1b1': ['kmin_mqml', 'kmax_mqml'], 
                'kqtl11.l1b1': ['kmin_mqtli','kmax_mqtli'], 
                'kqt12.l1b1': ['kmin_mqt','kmax_mqt'], 
                'kqt13.l1b1': ['kmin_mqt','kmax_mqt'],
                'kq4.l5b1': ['kmin_mqy_4.5k','kmax_mqy_4.5k'],
                'kq4.r5b1': ['kmin_mqy_4.5k','kmax_mqy_4.5k'],
                'kq5.l5b1': ['kmin_mqml_4.5k','kmax_mqml_4.5k'],
                'kq5.r5b1': ['kmin_mqml_4.5k','kmax_mqml_4.5k'],
                'kq6.l5b1': ['kmin_mqml_4.5k','kmax_mqml_4.5k'],
                'kq6.r5b1': ['kmin_mqml_4.5k','kmax_mqml_4.5k'],
                'kq7.l5b1': ['kmin_mqm_4.5k','kmax_mqm_4.5k'],
                'kq7.r5b1': ['kmin_mqm_4.5k','kmax_mqm_4.5k'],
                'kq8.l5b1': ['kmin_mqml','kmax_mqml'],
                'kq8.r5b1': ['kmin_mqml','kmax_mqml'],
                'kq9.l5b1': ['kmin_mqml','kmax_mqml'],
                'kq9.r5b1': ['kmin_mqml','kmax_mqml'],
                'kq10.l5b1': ['kmin_mqml', 'kmax_mqml'], 
                'kq10.r5b1': ['kmin_mqml', 'kmax_mqml'], 
                'kqtl11.l5b1': ['kmin_mqtli','kmax_mqtli'], 
                'kqt12.l5b1': ['kmin_mqt','kmax_mqt'], 
                'kqt13.l5b1': ['kmin_mqt','kmax_mqt'],
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
        if 'l1' in ii:
                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip1.b1']/350
        if 'l5' in ii:
                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip5.b1']/350


# %%
reset_delta_k(my_k_list)
line.vars['i_wire_ip1.b1'] = 350.0
line.vars['i_wire_ip5.b1'] = 0.0

tw_b1[['betx','alfx','bety','alfy'],'ip1']
tw_b1[['betx','alfx','bety','alfy'],'ip5']
tw_b1 = line.twiss(method='4d')

print('Before matching')
print(tw_b1.qx, tw_b1.qy)

variables_list = ['kq5.l1b1_delta', 'kq5.r1b1_delta']
opt = line.match(
    solve=False, # <- prepare the match without running it
    method='4d',
    vary=xt.VaryList(variables_list, limits=(-.001, 0.001),  step=1e-7),
    targets=xt.TargetList(qx=qx0, qy=qy0, tol=1e-7, tag='tunes'),
)

opt.step(10)

tw_b1 = line.twiss(method='4d')
print('After matching')
print(tw_b1.qx, tw_b1.qy)
# %%
plt.plot(tw_b1.s, (tw_b1.betx-betx_0)/betx_0, label='betx')
plt.plot(tw_b1.s, (tw_b1.bety-bety_0)/bety_0, label='betx')
# set xticks only at the IPs
plt.xticks([tw_b1['s','ip1'], tw_b1['s','ip2'], tw_b1['s','ip5'], tw_b1['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])
# %%

for ii in [1]:
        plt.figure()
        s0 = tw_b1['s',f'ip{ii}']

        plt.plot(tw_b1.s-s0, (tw_b1.betx-betx_0)/betx_0, '.-', label='betx')
        plt.plot(tw_b1.s-s0, (tw_b1.bety-bety_0)/bety_0, '.-', label='betx')
        plt.xlim(-5,5)
        plt.grid()
# %%
print('For 350 A and 8 sigma_y at TCTPV.4L1.B1\n')
print_k_summary(my_k_list)



# %%

reset_delta_k(my_k_list)
line.vars['i_wire_ip1.b1'] = 0.0
line.vars['i_wire_ip5.b1'] = 350.0

tw_b1[['betx','alfx','bety','alfy'],'ip1']
tw_b1[['betx','alfx','bety','alfy'],'ip5']
tw_b1 = line.twiss(method='4d')

print('Before matching')
print(tw_b1.qx, tw_b1.qy)


variables_list = ['kq4.l5b1_delta', 'kq4.r5b1_delta']
opt = line.match(
    solve=False, # <- prepare the match without running it
    method='4d',
    vary=xt.VaryList(variables_list, limits=(-.001, 0.001),  step=1e-7),
    targets=xt.TargetList(qx=qx0, qy=qy0, tol=1e-7, tag='tunes'),
)
opt.step(10)

tw_b1 = line.twiss(method='4d')
print('After matching')
print(tw_b1.qx, tw_b1.qy)

# %%
plt.plot(tw_b1.s, (tw_b1.betx-betx_0)/betx_0, label='$\Delta\\beta_x$-beating')
plt.plot(tw_b1.s, (tw_b1.bety-bety_0)/bety_0, label='$\Delta\\beta_x$-beating')
# set xticks only at the IPs
plt.legend()
plt.title('Wires at IP5@350 A, TCTH.4L5.B1@8$\sigma$, Q4 corrected')
plt.xticks([tw_b1['s','ip1'], tw_b1['s','ip2'], tw_b1['s','ip5'], tw_b1['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])
# %%
s0 = tw_b1['s','ip1']
plt.plot(tw_b1.s-s0, (tw_b1.betx-betx_0)/betx_0, '.-', label='$\Delta\\beta_x$-beating')
plt.plot(tw_b1.s-s0, (tw_b1.bety-bety_0)/bety_0, '.-', label='$\Delta\\beta_y$-beating')
plt.xlim(-5,5)
plt.grid()
# %%
print('For 350 A and 8 sigma_x at TCTPH.4L5.B1')
print_k_summary(my_k_list)

# %%
reset_delta_k(my_k_list)
line.vars['i_wire_ip1.b1'] = 350.0
line.vars['i_wire_ip5.b1'] = 0.0


tw_b1[['betx','alfx','bety','alfy'],'ip1']
tw_b1[['betx','alfx','bety','alfy'],'ip5']
tw_b1 = line.twiss(method='4d')

print('Before matching')
print(tw_b1.qx, tw_b1.qy)

variables_list = ['kq5.l1b1_delta', 'kq6.l1b1_delta']
opt = line.match(
    solve=False, # <- prepare the match without running it
    method='4d',
    vary=xt.VaryList(variables_list, limits=(-.001, 0.001),  step=1e-7),
    targets=xt.TargetList(qx=qx0, qy=qy0, tol=1e-7, tag='tunes'),
)


opt.step(10)

tw_b1 = line.twiss(method='4d')
print('After matching')
print(tw_b1.qx, tw_b1.qy)
# %%
plt.plot(tw_b1.s, (tw_b1.betx-betx_0)/betx_0, label='betx')
plt.plot(tw_b1.s, (tw_b1.bety-bety_0)/bety_0, label='betx')
# set xticks only at the IPs
plt.xticks([tw_b1['s','ip1'], tw_b1['s','ip2'], tw_b1['s','ip5'], tw_b1['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])

for ii in [1]:
        plt.figure()
        s0 = tw_b1['s',f'ip{ii}']

        plt.plot(tw_b1.s-s0, (tw_b1.betx-betx_0)/betx_0, '.-', label='betx')
        plt.plot(tw_b1.s-s0, (tw_b1.bety-bety_0)/bety_0, '.-', label='betx')
        plt.xlim(-5,5)
        plt.grid()

print('For 350 A and 8 sigma_y at TCTPV.4L1.B1')
print_k_summary(my_k_list)


# %% A test
reset_delta_k(my_k_list)
line.vars['i_wire_ip1.b1'] = 0.0
line.vars['i_wire_ip5.b1'] = 0.0
tw_ref = line.twiss(method='4d')

print('Reference')
print(tw_ref.qx, tw_ref.qy)



# %%
reset_delta_k(my_k_list)
line.vars['i_wire_ip1.b1'] = 350.0
line.vars['i_wire_ip5.b1'] = 0.0

tw_b1[['betx','alfx','bety','alfy'],'ip1']
tw_b1[['betx','alfx','bety','alfy'],'ip5']
tw_b1 = line.twiss(method='4d')

print('Before matching')
print(tw_b1.qx, tw_b1.qy)
variables_list = ['kq5.l1b1_delta', 
                  #'kq6.l1b1_delta',
                  'kq7.l1b1_delta',
                  #'kq8.l1b1_delta',
                  #'kq9.l1b1_delta',
                  #'kq10.l1b1_delta',
                  #'kqtl11.l1b1_delta',
                  #'kqt12.l1b1_delta',
                  #'kqt13.l1b1_delta',                  
                  ]

opt = line.match(
    solve=False, # <- prepare the match without running it
    start='s.ds.l1.b1', end='ip1',
    init=tw_ref, init_at=xt.START,
    method='4d',
    vary=
        xt.VaryList(variables_list, limits=(-.001, 0.001),  step=1e-7),
    targets=[
        xt.TargetSet(['mux', 'muy'], value=tw_ref, at=xt.END),
    ]
)
opt.step(10)

tw_b1 = line.twiss(method='4d')
print('After matching')
print(tw_b1.qx, tw_b1.qy)

plt.plot(tw_b1.s, (tw_b1.betx-tw_ref.betx)/tw_ref.betx, label='betx')
plt.plot(tw_b1.s, (tw_b1.bety-tw_ref.bety)/tw_ref.bety, label='betx')


# set xticks only at the IPs
plt.xticks([tw_b1['s','ip1'], tw_b1['s','ip2'], tw_b1['s','ip5'], tw_b1['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])

for ii in [1]:
        plt.figure()
        s0 = tw_b1['s',f'ip{ii}']
        plt.plot(tw_b1.s-s0, (tw_b1.betx-tw_ref.betx)/tw_ref.betx, '.-', label='betx')
        plt.plot(tw_b1.s-s0, (tw_b1.bety-tw_ref.bety)/tw_ref.bety, '.-', label='betx')
        plt.xlim(-700,5)
        plt.grid()
        plt.figure()
        plt.plot(tw_b1.s-s0, tw_b1.betx, label='betx')
        plt.plot(tw_b1.s-s0, tw_b1.bety, label='bety')
        plt.grid()
        plt.xlim(-700,5)
        plt.figure()
        plt.plot(tw_b1.s-s0, tw_b1.dx, '.-', label='dx')
        plt.plot(tw_b1.s-s0, tw_b1.dy, '.-', label='dy')
        plt.plot(tw_b1.s-s0, tw_ref.dx, 'sr', label='dx')
        plt.plot(tw_b1.s-s0, tw_ref.dy, 'sr', label='dy')
        plt.xlim(-700,5)


print('For 350 A and 8 sigma_y at TCTPV.4L1.B1')
print_k_summary(my_k_list)

# %%

reset_delta_k(my_k_list)

match_ip1 = True

if match_ip1:
        line.vars['i_wire_ip1.b1'] = 350.0
        line.vars['i_wire_ip5.b1'] = 0.0

        tw_b1 = line.twiss(method='4d')
        print('Before matching')
        print(tw_b1.qx, tw_b1.qy)

        variables_list = ['kq5.l1b1_delta', 
                        'kq6.l1b1_delta',
                        'kq7.l1b1_delta',
                        'kq8.l1b1_delta',
                        'kq9.l1b1_delta',
                        'kq10.l1b1_delta',
                        'kqtl11.l1b1_delta',
                        'kqt12.l1b1_delta',
                        #'kqt13.l1b1_delta',                  
                        ]
        opt = line.match(
        solve=False, # <- prepare the match without running it
        start='s.ds.l1.b1', end='ip1',
        init=tw_ref, init_at=xt.START,
        method='4d',
        vary=
                xt.VaryList(variables_list, limits=(-.001, 0.001),  step=1e-7),
        targets=[
                xt.TargetSet(['mux', 'muy','betx','bety','alfx','alfy','dx','dpx'], value=tw_ref, at=xt.END),
        ])       
        my_ip = 1
        my_tct = 'tctpv.4l1.b1'
        my_i = 'i_wire_ip1.b1'

else:
        line.vars['i_wire_ip1.b1'] = 0.0
        line.vars['i_wire_ip5.b1'] = 350.0
        tw_b1 = line.twiss(method='4d')
        print('Before matching')
        print(tw_b1.qx, tw_b1.qy)
        variables_list = ['kq4.l5b1_delta', 
                        'kq5.l5b1_delta', 
                        #'kq6.l5b1_delta',
                        'kq7.l5b1_delta',
                        'kq8.l5b1_delta',
                        'kq9.l5b1_delta',
                        'kq10.l5b1_delta',
                        'kqtl11.l5b1_delta',
                        'kqt12.l5b1_delta',
                        #'kqt13.l5b1_delta',                  
                        ]
        opt = line.match(
        solve=False, # <- prepare the match without running it
        start='s.ds.l5.b1', end='ip5',
        init=tw_ref, init_at=xt.START,
        method='4d',
        vary=
                xt.VaryList(variables_list, limits=(-.001, 0.001),  step=1e-7),
        targets=[
                xt.TargetSet(['mux', 'muy','betx','bety','alfx','alfy','dx','dpx'], value=tw_ref, at=xt.END),
        ])
        my_ip = 5
        my_tct = 'tctph.4l5.b1'
        my_i = 'i_wire_ip5.b1'

opt.step(10)

tw_b1 = line.twiss(method='4d')
print('After matching')
print(tw_b1.qx, tw_b1.qy)

plt.plot(tw_b1.s, (tw_b1.betx-tw_ref.betx)/tw_ref.betx, label='$\Delta\\beta_x/\\beta_{x0}$')
plt.plot(tw_b1.s, (tw_b1.bety-tw_ref.bety)/tw_ref.bety, label='$\Delta\\beta_y\\beta_{y0}$')

# set xticks only at the IPs
plt.xticks([tw_b1['s','ip1'], tw_b1['s','ip2'], tw_b1['s','ip5'], tw_b1['s','ip8']], ['IP1', 'IP2', 'IP5', 'IP8'])
plt.xlabel('s along the ring [no units]')
plt.ylabel('relative $\\beta$-beating')
plt.title(f'Wire at IP{my_ip}@350 A, {my_tct}@8$\sigma@30cm$ with optics at {my_optics} cm')
plt.grid()
plt.legend(loc='lower right')

for ii in [my_ip]:
        plt.figure()
        s0 = tw_b1['s',f'ip{ii}']
        plt.plot(tw_b1.s-s0, (tw_b1.betx-tw_ref.betx)/tw_ref.betx, '.-', label='$\Delta\\beta_x/\\beta_{x0}$')
        plt.plot(tw_b1.s-s0, (tw_b1.bety-tw_ref.bety)/tw_ref.bety, '.-', label='$\Delta\\beta_y\\beta_{y0}$')
        plt.xlim(-700,5)
        plt.grid()
        plt.xlabel(f's from IP{my_ip} [m]')
        plt.ylabel('relative $\\beta$-beating')
        plt.title(f'Wire at IP{my_ip}@350 A, {my_tct}@8$\sigma@30cm$ with optics at {my_optics} cm')
        plt.legend(loc='upper left')

        

        # set the elements_list taking the list line.element_names from 's.ds.l1.b1' to 'ip1'
        element_list = line.element_names[line.element_names.index('s.ds.l1.b1'):line.element_names.index('ip1')+1]

        s_ref = line.get_s_position('ip1')

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
        plt.plot(tw_b1.s-s0, tw_b1.betx, 'r', label='$\\beta_x$')
        plt.plot(tw_b1.s-s0, tw_b1.bety, 'b', label='$\\beta_y$')
        plt.ylabel('$\\beta_x$ and $\\beta_y$ [m]')
        plt.legend(loc='lower right')

        plt.figure()
        plt.plot(tw_b1.s-s0, tw_b1.dx, '.-b', label='dx')
        plt.plot(tw_b1.s-s0, tw_b1.dy, '.-r', label='dy')
        plt.plot(tw_b1.s-s0, tw_ref.dx, 'ob', label='dx$_0$')
        plt.plot(tw_b1.s-s0, tw_ref.dy, 'or', label='dy$_0$')
        plt.xlim(-700,5)
        plt.ylim(-0.5, 2.5)
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
        line.vars['i_wire_ip1.b1'] = my_current
        line.vars['i_wire_ip5.b1'] = 0.0
else:
        line.vars['i_wire_ip1.b1'] = 0.0
        line.vars['i_wire_ip5.b1'] = my_current

tw_b1 = line.twiss(method='4d')

print(tw_b1.qx-tw_ref.qx, tw_b1.qy-tw_ref.qy)
for ii in [my_ip]:
        plt.figure()
        s0 = tw_b1['s',f'ip{ii}']
        plt.plot(tw_b1.s-s0, (tw_b1.betx-tw_ref.betx)/tw_ref.betx, '.-', label='$\Delta\\beta_x/\\beta_{x0}$')
        plt.plot(tw_b1.s-s0, (tw_b1.bety-tw_ref.bety)/tw_ref.bety, '.-', label='$\Delta\\beta_y\\beta_{y0}$')
        #plt.xlim(-700,5)
        plt.grid()
        plt.xlabel(f's from IP{my_ip} [m]')
        plt.ylabel('relative $\\beta$-beating')
        plt.title(f'Wire at IP{my_ip}@{my_current} A, {my_tct}@8$\sigma@30cm$ with optics at {my_optics} cm')
        plt.legend(loc='upper left')
        plt.ylim(-.028, 0.028)
# %%    
        

        # # set the elements_list taking the list line.element_names from 's.ds.l1.b1' to 'ip1'
        # element_list = line.element_names[line.element_names.index('s.ds.l1.b1'):line.element_names.index('ip1')+1]

        # s_ref = line.get_s_position('ip1')

        # for element in element_list:
        #         # if line.element_dict[element] belongs to xt.Quadrupole
        #         if line[element].__class__.__name__ == 'Quadrupole':
        #                 my_dict = line[element].to_dict()
        #                 my_s = line.get_s_position(element)-s_ref
        #                 my_length = my_dict['length']
        #                 my_k = my_dict['k1']
        #                 # plot a red filled rectangle starting at my_s, 0 with length my_length and height my_k 
        #                 plt.gca().add_patch(plt.Rectangle((my_s, 0), my_length, my_k, color='r', alpha=0.5))
                       


# gain = 1
# if my_ip == 1:
#         line.vars['i_wire_ip1.b1'] = 350.0
#         line.vars['i_wire_ip5.b1'] = 0.0
# tw_100_on = line.twiss(method='4d')

# list_betx_beating = []
# list_bety_beating = []
# list_gain = np.linspace(0, 1, 11)
# for gain in list_gain:
#         print(gain)
#         for ii in variables_list:
#                 line.vars[ii] = line.vars[ii] * gain
#         if my_ip == 1:
#                 line.vars['i_wire_ip1.b1'] = 350.0 * gain
#                 line.vars['i_wire_ip5.b1'] = 0.0 * gain
#         tw_aux = line.twiss(method='4d')

#         list_betx_beating.append(np.max(np.abs(tw_aux.betx-tw_100_on.betx)/tw_100_on.betx))
#         list_bety_beating.append(np.max(np.abs(tw_aux.bety-tw_100_on.bety)/tw_100_on.bety))

        

# %%
