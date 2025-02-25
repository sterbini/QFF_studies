# %%
import json
from matplotlib import pyplot as plt
import xtrack as xt
import numpy as np

my_beam = 'b1' # this is the beam 4 in madx sense

for my_optics in [22, 24,26, 28, 30]:
       print(8*'*', f'Optics: {my_optics} cm', 8*'*')
       for my_ip in [1, 5]:
                if my_ip==1:
                        x_position = -15000
                        my_file='knob_dict_350A_8sigma@30cm_ip1_beta30_b1.json'
                else:
                        x_position = 10000
                        my_file='knob_dict_350A_8sigma@30cm_ip5_beta30_b1.json'

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

                line = collider[f'lhc{my_beam}']

                wire_s_ip1 = line.get_s_position('tctpv.4l1.b1_entry')-(line.get_s_position('tctpv.4l1.b1_entry')-line.get_s_position('tctpv.4l1.b1_exit'))/2
                wire_s_ip5 = line.get_s_position('tctph.4l5.b1_entry')-(line.get_s_position('tctph.4l5.b1_entry')-line.get_s_position('tctph.4l5.b1_exit'))/2 

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
                                at_s=wire_s_ip1)
                # bottom wire in IR1
                line.insert_element(name=f'bbwc.b.4l1.{my_beam}',
                                element=xt.Wire(
                                        L_phy=1, 
                                        L_int=2, 
                                        current=0.0, 
                                        xma=0.0, 
                                        yma=-1  # very far from the beam
                                        ),
                                at_s=wire_s_ip1)

                # upper wire in IR5
                line.insert_element(name=f'bbwc.e.4l5.{my_beam}',
                                element=xt.Wire(
                                        L_phy=1, 
                                        L_int=2,
                                        current= 0.0,
                                        xma=1, # very far from the beam
                                        yma=0.0 
                                        ),
                                at_s=wire_s_ip5)
                # bottom wire in IR5
                line.insert_element(name=f'bbwc.i.4l5.{my_beam}',
                                element=xt.Wire(
                                        L_phy=1, 
                                        L_int=2, 
                                        current=0.0, 
                                        xma=-1, # very far from the beam
                                        yma=0.0  
                                        ),
                                at_s=wire_s_ip5)

                s_ip = line.get_s_position(at_elements=f'ip{my_ip}')

                line.vars['co_y_wire_ip1.b1'] = 0
                line.vars['co_x_wire_ip1.b1'] = 0
                line.vars['co_y_wire_ip5.b1'] = 0
                line.vars['co_x_wire_ip5.b1'] = 0        

                line.element_refs['bbwc.t.4l1.b1'].current = line.vars['i_wire_ip1.b1']
                line.element_refs['bbwc.t.4l1.b1'].yma = line.vars['d_wire_ip1.b1'] + line.vars['co_y_wire_ip1.b1']

                line.element_refs['bbwc.b.4l1.b1'].current = line.vars['i_wire_ip1.b1'] 
                line.element_refs['bbwc.b.4l1.b1'].yma = -line.vars['d_wire_ip1.b1'] + line.vars['co_y_wire_ip1.b1']

                line.element_refs['bbwc.e.4l5.b1'].current = line.vars['i_wire_ip5.b1']
                line.element_refs['bbwc.e.4l5.b1'].xma = line.vars['d_wire_ip5.b1'] + line.vars['co_x_wire_ip5.b1']

                line.element_refs['bbwc.i.4l5.b1'].current = line.vars['i_wire_ip5.b1']
                line.element_refs['bbwc.i.4l5.b1'].xma = -line.vars['d_wire_ip5.b1'] + line.vars['co_x_wire_ip5.b1']

                for ii in [2,8]:
                        collider.vars[f'on_x{ii}h'] = 0.0
                        collider.vars[f'on_x{ii}v'] = 0.0
                        collider.vars[f'on_sep{ii}h'] = 0.0
                        collider.vars[f'on_sep{ii}v'] = 0.0

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

                # import json file 
                with open(my_file) as f:
                        data = json.load(f)

                line.vars['d_wire_ip1.b1'] = data['tct_opening_in_sigma'] * data['sigma_y_at_tctpv_4l1_b1'] + data['wire_retraction']
                line.vars['d_wire_ip5.b1'] = data['tct_opening_in_sigma'] * data['sigma_x_at_tctph_4l5_b1'] + data['wire_retraction']



                for ii in data['k_0']:
                        assert data['k_0'][ii] == line.vars[ii]._get_value()
                
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

                def reset_delta_k(my_k_list):
                        for kk in my_k_list:
                                collider.vars[f'{kk}_delta'] = 0.000000

                reset_delta_k(my_k_list)

                for ii in my_k_list:
                        collider.vars[f'{ii}_0'] = collider.vars[ii]._get_value()
                        collider.vars[f'{ii}_delta'] = 0.000000
                        if 'r1' in ii:
                                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip1.b1']/350
                        if 'l1' in ii:
                                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip1.b1']/350
                
                        if 'r5' in ii:
                                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip5.b1']/350
                        if 'l5' in ii:
                                collider.vars[ii] = collider.vars[f'{ii}_0'] + collider.vars[f'{ii}_delta']*collider.vars['i_wire_ip5.b1']/350


                for ii in data['k_delta']:
                        collider.vars[f'{ii}_delta'] = data['k_delta'][ii]

                tw_ref = line.twiss(method='4d')


                line.vars['i_wire_ip1.b1'] = data['i_wire_ip1.b1']
                line.vars['i_wire_ip5.b1'] = data['i_wire_ip5.b1']

                tw_new = line.twiss(method='4d')
                s0 = tw_new['s',f'ip{my_ip}']
                plt.figure()
                plt.plot(tw_ref['s']-s0, (tw_new['betx']- tw_ref['betx'])/tw_ref['betx'], label='$\Delta\\beta_x/\\beta_{x0}$')
                plt.plot(tw_ref['s']-s0, (tw_new['bety']- tw_ref['bety'])/tw_ref['bety'], label='$\Delta\\beta_y/\\beta_{y0}$')

                t = plt.text(x_position, -0.025, f'''$\Delta Q_x$ = {tw_new.qx-tw_ref.qx:.2e}
$\Delta Q_y$ = {tw_new.qy-tw_ref.qy:.2e}''', fontsize=12)
                t.set_bbox(dict(facecolor='white', edgecolor='red'))

                plt.grid()
                plt.xlabel(f's from IP{my_ip} [m]')
                plt.ylabel('relative $\\beta$-beating')
                plt.title(f'Wire at IP{my_ip}@350 A, TCT@8$\sigma$ with optics at $\\beta$*={my_optics} cm')
                plt.legend(loc='upper left')
                plt.ylim(-0.055, 0.055)
                plt.savefig(f'plots/knob_{my_beam}_ip{my_ip}_{my_optics}cm.png')

# %%
