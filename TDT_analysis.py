# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 16:01:52 2022

@author: Dkeum
"""

import tdt
import os
import keyboard as kb
import matplotlib.pyplot as plt  
import numpy as np  
import pandas as pd
import openpyxl
import peakdetect as pdetec
import scipy.stats as stats





def VCU_spike_analysis(tank, block, sort_name = 'autosort1'):
    trough_thres = 0.2e-3 #s
    peak_thres = 1.0e-3 #s
    #Read a block
    VCU_data_path = 'D:/OneDrive - University of Maryland School of Medicine/in_vivo_data'
    source = '/' + tank + '/' + 'Block-'+ str(block)
    path = VCU_data_path + source
    data = tdt.read_block(path, sortname = sort_name)
    fs = data.snips.Snip.fs
    
    #
    plt.ioff()
    
    #Convert TDT data to DataFrame
    x_point = [t for t in range(0, 34)]
    x_time = np.array(x_point) * 1/fs
    col = ['ts', 'Channel', 'Sort'] + x_point
    d = np.column_stack([data.snips.Snip.ts, data.snips.Snip.chan, data.snips.Snip.sortcode, data.snips.Snip.data])
    
    df = pd.DataFrame(data = d, columns = col)
    #df.to_csv('VCU.csv')
    
    drop_index  = df.index[(df['Channel'] > 32) | (df['Channel'] < 0) | (df['Sort'] > 30)]
    df.drop(index = drop_index, inplace = True)
    df.reset_index(inplace = True)
    #Add CellID
    df = df.astype({'Channel': int, 'Sort': int})
    df["CellID"] = df['Channel'].astype(str) + '_' + df['Sort'].astype(str)
    
    
    
    #Buffer for storing avg spike of each cell
    col = ['CellID'] + x_point
    df_AVG_spike = pd.DataFrame(columns = col)
    cells =  pd.unique(df['CellID'])
    cells = np.sort(cells)
    #Buffer for storing average spike width of each cell
    col = ['Tank', 'Block', 'CellID', 'P2P_width', 'P1_T_width', 'P1_T_ratio', 'P2_T_width', 'P2_T_ratio']
    df_spike_property = pd.DataFrame(columns = col)
    df_spike_property_AVG = pd.DataFrame(columns = col)
    biphasic_spikes = pd.DataFrame(columns = ['Tank', 'Block', 'CellID', 'Index'] + x_point)
    biphasic_spike_avg = pd.DataFrame(columns = ['Tank', 'Block', 'CellID'] + x_point)
        
        
    for this_cell in cells:
        this_cell_spikes = df[df['CellID'] == this_cell][x_point]

        
        #Temporary code for displaying spike shapes###############
        #Find and gather bipolar spikes 
        #print(this_cell)
        width = np.nan
        isBiphasic = False
        cell_spike_width = np.array([])
        
        print(this_cell)
        for index, row in this_cell_spikes.iterrows():
            y_sig = this_cell_spikes.loc[index].to_numpy()
            extrema = pdetec.peakdetect(y_sig, x_time, 10)
            peaks, troughs = np.array(extrema, dtype = object)
            #Check the spike biphasic. if so, obtain width between trough and peak (Merrikhi et al., 2021)
            criteria1 = False
            criteria2 = False
            criteria3 = False
            
            if len(troughs) > 0 and len(peaks) > 1:
                #first basin between 0.0 and 0.2ms
                idx_trough = 0
                idx_peak = 0
                
                #Criteria 1: Is the first peak is no sooner than 0.2 ms ?
                if peaks[0][0] >= trough_thres:
                    criteria1 = True
                    idx_trough = 0

                #Criteria 2: Is the 2nd peak is no later than 1.0 ms ?
                if peaks[1][0] <= peak_thres:
                    criteria2 = True
                        
                #Criteria 3: Is the 1st trough is between 2 peaks?
                if troughs[0][0] > peaks[0][0] and troughs[0][0] < peaks[1][0]:
                    criteria3 = True
                
                if criteria1 & criteria2 & criteria3:
                    isBiphasic = True
                    #spike_width = peaks[idx_peak][0] - troughs[idx_trough][0]                    
                    peak_to_peak_width = peaks[1][0] - peaks[0][0]
                    peak1_trough_width = troughs[0][0] - peaks[0][0]
                    peak1_trough_ratio = peaks[0][1] / troughs[0][1]
                    peak2_trough_width = peaks[1][0] - troughs[0][0]
                    peak2_trough_ratio = peaks[1][1] / troughs[0][1]
                    
                    biphasic_spikes.loc[len(biphasic_spikes)] = [tank, block, this_cell, index] + list(y_sig)
                    df_spike_property.loc[len(df_spike_property)] = [tank, block, this_cell, \
                                                                     peak_to_peak_width, peak1_trough_width, peak1_trough_ratio,\
                                                                         peak2_trough_width, peak2_trough_ratio]
                    #print(spike_width)
                
        if isBiphasic:    
            #width = np.mean(cell_spike_width)
            this_cell_bi_spikes = biphasic_spikes[biphasic_spikes['CellID'] == this_cell][x_point]
            AVG_biphasic_spike = this_cell_bi_spikes.mean(axis = 0).to_list()
            biphasic_spike_avg.loc[len(biphasic_spike_avg)] = [tank, block, this_cell] + AVG_biphasic_spike
            this_cell_bi_spikes = biphasic_spikes[biphasic_spikes['CellID'] == this_cell][x_point]
            for bi_index, bi_row in this_cell_bi_spikes.iterrows():
                plt.plot(x_time, this_cell_bi_spikes.loc[bi_index], color = 'gray')
            #print(this_cell_bi_spikes)
            plt.plot(x_time, AVG_biphasic_spike, linewidth = 2.0, color = 'red')
            plt.xticks([0,3e-4,6e-4,9e-4,12e-4])
            plt.ticklabel_format(axis = 'x', style = 'scientific', scilimits = (-3, -3), useMathText = False)
            figtitle = filename = str(tank) +"_"+ str(block) +"_"+ this_cell
            plt.title(figtitle)
            plt.savefig('test_fig/' + filename)
            plt.cla()
            #Averaging spike property
            this_cell_spike_property = df_spike_property[df_spike_property['CellID'] == this_cell]
            #avg_properties = this_cell_spike_property.iloc[:,3:8].mean(axis = 0)
            #print(avg_properties)
            avg_properties = this_cell_spike_property.iloc[:,3:8].mean(axis = 0).to_list()
            df_spike_property_AVG.loc[len(df_spike_property_AVG)] = [tank, block, this_cell] + avg_properties
            
            
        #df_spike_width.loc[len(df_spike_width)] = [tank, block, this_cell, width]
        
        #print(this_cell + '_' + str(width))
            
            
        '''
            plt.plot(x_time, y_sig)
            print(peak)
            #plt.plot(peak)
        
        
        plt.xticks([0,3e-4,6e-4,9e-4,12e-4])
        plt.ticklabel_format(axis = 'x', style = 'scientific', scilimits = (-3, -3), useMathText = False)
        figtitle = filename = str(tank) +"_"+ str(block) +"_"+ this_cell
        plt.title(figtitle)
        filename = figtitle+'.png'
        plt.savefig('test_fig/' + filename)
        plt.cla()
        '''
        #########################################################
        #print('Press n for next cell.')
        #kb.wait('n')
        AVG_spike = this_cell_spikes.mean(axis = 0).to_list()
        df_AVG_spike.loc[len(df_AVG_spike)] = [this_cell] + AVG_spike
        
        
        
        
    #prinret(df_AVG_spike)
    return(biphasic_spike_avg, df_spike_property_AVG)
    
    
def Burst_detection():
    pass



def Proc_MS_freq(tank, block, sort_name = 'test', figures = False):
    #Read a block
    source = './' + tank + '/' + 'Block-'+ str(block)
    data = tdt.read_block(source, sortname = sort_name)
    fs = data.snips.eNeu.fs
    
    #Get stimuli time 
    som_mask = data.epocs.ssom.data.astype(bool)
    vis_mask = data.epocs.svis.data.astype(bool)
    ts = data.epocs.ssom.onset
    
    t_som = ts[som_mask]
    t_vis = ts[vis_mask]
    t_bim = np.intersect1d(t_som, t_vis)
    t_som = np.setdiff1d(t_som, t_bim)
    t_vis = np.setdiff1d(t_vis, t_bim)

    t_vis = t_vis + 0.1
    t_som = t_som + 0.25
    t_bim = t_bim + 0.25
    
    #Making data frame for analysis
    
    ts_col = ['Tank', 'Block', 'CellID', 'Cell_spon_ts', 'Cell_resp_ts', 'Vis_spon_ts', \
              'Vis_resp_ts', 'Som_spon_ts', 'Som_resp_ts', 'Bim_spon_ts', 'Bim_resp_ts']
    Block_result_ts = pd.DataFrame(columns = ts_col)
    sp_col = ['Tank', 'Block', 'CellID', 'Cell_spon_avg', 'Cell_resp_avg','Vis_spon_avg', \
              'Vis_resp_avg', 'Som_spon_avg', 'Som_resp_avg', 'Bim_spon_avg', 'Bim_resp_avg']
    Block_result_sp = pd.DataFrame(columns = sp_col)
    ft_col = ['Tank', 'Block', 'CellID','Cell_spon__ISI_median', 'Cell_resp_ISI_median']
    Block_result_ft = pd.DataFrame(columns = ft_col)
    
    CellID = []
    drop_ind = []
    
    for ind in range(len(data.snips.eNeu.sortcode)):
        id = str(data.snips.eNeu.chan[:,0][ind])+'_'+\
            str(data.snips.eNeu.sortcode[:,0][ind]) 
        CellID.append(id)
        if data.snips.eNeu.chan[:,0][ind] < 0 or data.snips.eNeu.chan[:,0][ind] > 16:
            drop_ind.append(ind)
        if data.snips.eNeu.sortcode[:,0][ind] > 30 or data.snips.eNeu.sortcode[:,0][ind] < 1:
            drop_ind.append(ind)

        
    CellID = np.array(CellID)
    spike_time = [str(a) for a in range(30)] 
    d_col = ['Time', 'CellID']
    d = np.column_stack([data.snips.eNeu.ts, CellID, data.snips.eNeu.data])
    NeuAct = pd.DataFrame(data = d, columns = d_col+spike_time)
    
    NeuAct = NeuAct.astype({'Time': 'float'})
    #Drop cells with error
    drop_ind = list(set(drop_ind))
    NeuAct = NeuAct.drop(drop_ind)
    NeuAct.reset_index(drop = True)
    cells = pd.unique(NeuAct['CellID'])
    
    
    for cell in cells:
        cell_spontaneous_ts = np.array([])
        cell_response_ts = np.array([])
        vis_spontaneous_ts = np.array([])
        som_spontaneous_ts = np.array([])
        bim_spontaneous_ts = np.array([])
        vis_response_ts = np.array([])
        som_response_ts = np.array([])
        bim_response_ts = np.array([])
        
        cell_spontaneous_spikes = pd.DataFrame()
        cell_response_spikes = pd.DataFrame()
        vis_spontaneous_spikes = pd.DataFrame()
        som_spontaneous_spikes = pd.DataFrame()
        bim_spontaneous_spikes = pd.DataFrame()
        vis_response_spikes = pd.DataFrame()
        som_response_spikes = pd.DataFrame()
        bim_response_spikes = pd.DataFrame()
        
        cell_spontaneous_spike_avg = np.array([])
        cell_response_spike_avg = np.array([])
        vis_spontaneous_spike_avg = np.array([])
        som_spontaneous_spike_avg = np.array([])
        bim_spontaneous_spike_avg = np.array([])
        vis_response_spike_avg = np.array([])
        som_response_spike_avg = np.array([])
        bim_response_spike_avg = np.array([])
        
        this_cell = NeuAct[NeuAct['CellID']==cell]
        
        #Visual response
        for index in range(len(t_vis)):
            stim = t_vis[index]
            pre = this_cell[(this_cell['Time'] >= stim - 1) & (this_cell['Time'] < stim )]
            ongoing = this_cell[(this_cell['Time'] >= stim) & (this_cell['Time'] < stim + 1)]
            spontaneous = pre['Time'].to_numpy()
            response = ongoing['Time'].to_numpy()
            
            if len(response) > 0:
                if len(vis_response_ts) == 0:
                    response = response - response[0]
                    vis_response_ts = np.append(vis_response_ts, response)
                else:
                    response = response - response[0] + vis_response_ts[-1]
                    vis_response_ts = np.append(vis_response_ts, response[1:])
            if len(spontaneous) > 0:
                if len(vis_spontaneous_ts) == 0:
                    spontaneous = spontaneous - spontaneous[0]
                    vis_spontaneous_ts = np.append(vis_spontaneous_ts, spontaneous)
                    
                else:
                    spontaneous = spontaneous - spontaneous[0] + vis_spontaneous_ts[-1]
                    vis_spontaneous_ts = np.append(vis_spontaneous_ts, spontaneous[1:])
                    
            
            spon_spikes = pre.iloc[:, 2:32].astype('float')
            respon_spikes = ongoing.iloc[:, 2:32].astype('float')
            
            vis_spontaneous_spikes = pd.concat([vis_spontaneous_spikes, spon_spikes])
            vis_response_spikes = pd.concat([vis_response_spikes, respon_spikes])
            
        vis_spontaneous_spike_avg = vis_spontaneous_spikes.mean(axis = 0)
        vis_response_spike_avg = vis_response_spikes.mean(axis = 0)
        
        #Tactile (somatosensory) response
        for index in range(len(t_som)):
            stim = t_som[index]
            pre = this_cell[(this_cell['Time'] >= stim - 1) & (this_cell['Time'] < stim )]
            ongoing = this_cell[(this_cell['Time'] >= stim) & (this_cell['Time'] < stim + 1)]
            spontaneous = pre['Time'].to_numpy()
            response = ongoing['Time'].to_numpy()
            if len(response) > 0:
                if len(som_response_ts) == 0:
                    response = response - response[0]
                    som_response_ts = np.append(som_response_ts, response)
                else:
                    response = response - response[0] + som_response_ts[-1]
                    som_response_ts = np.append(som_response_ts, response[1:])
            if len(spontaneous) > 0:
                if len(som_spontaneous_ts) == 0:
                    spontaneous = spontaneous - spontaneous[0]
                    som_spontaneous_ts = np.append(som_spontaneous_ts, spontaneous)
                else:
                    spontaneous = spontaneous - spontaneous[0] + som_spontaneous_ts[-1]
                    som_spontaneous_ts = np.append(som_spontaneous_ts, spontaneous[1:])
                    
            spon_spikes = pre.iloc[:, 2:32].astype('float')
            respon_spikes = ongoing.iloc[:, 2:32].astype('float')
          
            som_spontaneous_spikes = pd.concat([som_spontaneous_spikes, spon_spikes])
            som_response_spikes = pd.concat([som_response_spikes, respon_spikes])
      
        som_spontaneous_spike_avg = som_spontaneous_spikes.mean(axis = 0)
        som_response_spike_avg = som_response_spikes.mean(axis = 0)

        #Bimodal response
        for index in range(len(t_bim)):
            stim = t_bim[index]
            pre = this_cell[(this_cell['Time'] >= stim - 1) & (this_cell['Time'] < stim )]
            ongoing = this_cell[(this_cell['Time'] >= stim) & (this_cell['Time'] < stim + 1)]
            spontaneous = pre['Time'].to_numpy()
            response = ongoing['Time'].to_numpy()
            if len(response) > 0:
                if len(bim_response_ts) == 0:
                    response = response - response[0]
                    bim_response_ts = np.append(bim_response_ts, response)
                else:
                    response = response - response[0] + bim_response_ts[-1]
                    bim_response_ts = np.append(bim_response_ts, response[1:])
            if len(spontaneous) > 0:
                if len(bim_spontaneous_ts) == 0:
                    spontaneous = spontaneous - spontaneous[0]
                    bim_spontaneous_ts = np.append(bim_spontaneous_ts, spontaneous)
                else:
                    spontaneous = spontaneous - spontaneous[0] + bim_spontaneous_ts[-1]
                    bim_spontaneous_ts = np.append(bim_spontaneous_ts, spontaneous[1:])
            
            spon_spikes = pre.iloc[:, 2:32].astype('float')
            respon_spikes = ongoing.iloc[:, 2:32].astype('float')
            
            bim_spontaneous_spikes = pd.concat([bim_spontaneous_spikes, spon_spikes])
            bim_response_spikes = pd.concat([bim_response_spikes, respon_spikes])
            
        bim_spontaneous_spike_avg = bim_spontaneous_spikes.mean(axis = 0)
        bim_response_spike_avg = bim_response_spikes.mean(axis = 0)
        
        
        cell_spon_length = 0
        if len(vis_spontaneous_ts)>0:
            cell_spontaneous_ts = vis_spontaneous_ts
            if len(som_spontaneous_ts) > 0:
                temp_ts = som_spontaneous_ts - som_spontaneous_ts[0] + cell_spontaneous_ts[-1]
                cell_spontaneous_ts = np.append(cell_spontaneous_ts, temp_ts[1:])
            if len(bim_spontaneous_ts) > 0:
                temp_ts = bim_spontaneous_ts - bim_spontaneous_ts[0] + cell_spontaneous_ts[-1]
                cell_spontaneous_ts = np.append(cell_spontaneous_ts, temp_ts[1:])
            
            cell_spon_length = cell_spontaneous_ts[-1]-cell_spontaneous_ts[0]
        
        cell_resp_length = 0
        if len(vis_response_ts)>0:
            cell_response_ts = vis_response_ts
            if(len(som_response_ts)>0):
                temp_ts = som_response_ts - som_response_ts[0] + cell_response_ts[-1]
                cell_response_ts = np.append(cell_response_ts, temp_ts[1:])
            if(len(bim_response_ts)>0):
                temp_ts = bim_response_ts - bim_response_ts[0] + cell_response_ts[-1]
                cell_response_ts = np.append(cell_response_ts, temp_ts[1:])
            
            cell_resp_length = cell_response_ts[-1]-cell_response_ts[0]
        
        '''
        print(cell)
        #print(len(cell_response_ts)/len(cell_spontaneous_ts))
        print('spontaneous')
        print(cell_spontaneous_ts)
        print(len(cell_spontaneous_ts))
        print('response')
        print(cell_response_ts)
        print(len(cell_response_ts))
        
        ts_col = ['Tank', 'Block', 'CellID', 'Cell_spon_ts', 'Cell_resp_ts', 'Vis_spon_ts', \
              'Vis_resp_ts', 'Som_spon_ts', 'Som_resp_ts', 'Bim_spon_ts', 'Bim_resp_ts']
                
                
        '''
        cell_spontaneous_spikes = pd.concat([vis_spontaneous_spikes, som_spontaneous_spikes, bim_spontaneous_spikes])
        cell_response_spikes = pd.concat([vis_response_spikes, som_response_spikes, bim_response_spikes])
        
        cell_spontaneous_spike_avg = cell_spontaneous_spikes.mean(axis = 0)
        cell_response_spike_avg = cell_response_spikes.mean(axis = 0)
        
        cell_spon_ISI = 0
        cell_resp_ISI = 0
        cell_spon_ISI_median = 0
        cell_resp_ISI_median = 0
        
        if(cell_spon_length > 1):
            cell_spon_ISI = cell_spontaneous_ts[1:] - cell_spontaneous_ts[:-1]
            cell_spon_ISI_median = np.median(cell_spon_ISI)
        if(cell_resp_length > 1):
            cell_resp_ISI = cell_response_ts[1:] - cell_response_ts[:-1]
            cell_resp_ISI_median = np.median(cell_resp_ISI)
            
        #Figure
        
        if (figures == True):
            plt.ioff()
            interval = 1/fs * 1000 #ms
            spike_x = np.array([a*interval for a in range(0, 30)])
            figtitle = tank + '_' + str(block) + '_' + cell
            fig, ax = plt.subplots(1, 2, figsize = (6,3))
            ax[0].plot(spike_x, cell_spontaneous_spike_avg, '-r')
            ax[0].plot(spike_x, cell_response_spike_avg, '-k')
            ax[0].set_xticks([0, 0.3, 0.6, 0.9, 1.2])
            ax[1].hist(cell_resp_ISI, bins = 100, range = [0E-3, 25E-3], color = 'black', label = 'Response')
            ax[1].hist(cell_spon_ISI, bins = 100, range = [0E-3, 25E-3], color = 'red', label = 'Spontaneous')
            ax[1].set_xticks([0,5e-3,10e-3,15e-3,20e-3,25e-3] )
            
            ax[1].ticklabel_format(axis = 'x', style = 'scientific', scilimits = (-3, -3), useMathText = False)
            fig.suptitle(figtitle)
            ax[1].legend(loc = 'upper right')
            filename = figtitle +'.png'
            plt.savefig('result_fig/' + filename)
                        
            
        #Container for result output
        result_ts = [tank, block, cell, cell_spontaneous_ts, cell_response_ts, \
                     vis_spontaneous_ts, vis_response_ts, som_spontaneous_ts, som_response_ts,\
                         bim_spontaneous_ts, bim_response_ts]
        result_sp = [tank, block, cell, cell_spontaneous_spike_avg, cell_response_spike_avg, \
                     vis_spontaneous_spike_avg, vis_response_spike_avg, som_spontaneous_spike_avg, som_response_spike_avg,\
                         bim_spontaneous_spike_avg, bim_response_spike_avg]
        result_ft = [tank, block, cell, cell_spon_ISI_median, cell_resp_ISI_median]
        Block_result_ts.loc[len(Block_result_ts)] = result_ts
        Block_result_sp.loc[len(Block_result_sp)] = result_sp
        Block_result_ft.loc[len(Block_result_ft)] = result_ft
        
    
    return [Block_result_ts, Block_result_sp, Block_result_ft]
                                             

def MS_freq_analysis(DataSet = 'MMP_Log.xlsx', Sheet = 'Final_block', Show_fig = True):
    FileList = pd.read_excel(DataSet, Sheet)
    Groups = pd.unique(FileList['Group'])
    #print(FileList)
    df_feature = pd.DataFrame()
    for grp in Groups:  
        this_group = FileList[FileList['Group'] == grp] 
        Treatments = pd.unique(this_group['Treatment'])
        for treat in Treatments: 
            this_treat = this_group[this_group['Treatment'] == treat]
            Animals = pd.unique(this_treat['Animal'])
            df_header = pd.DataFrame([[grp, treat]], columns = ['Group', 'Treatment'])
            for ani in Animals:
                this_animal = this_treat[this_treat['Animal'] == ani]
                Blocks = this_animal['Block'].to_numpy()
                
                #print(grp, treat, ani, Blocks )
                for block in Blocks:
                    #Get sort name 
                    sort_name = 'None'
                    try:
                        path = './' + ani + '/Block-' + str(block) + '/sort'
                        dirs = os.listdir(path)
                        for dir in dirs:
                            if os.path.isdir(path + '/' + dir):
                                sort_name = dir
                                break
                    except:
                        print('This block has no Bayesian sort..')
                    
                    print(grp, treat, ani, block, sort_name)
                    
                    try:
                        result = Proc_MS_freq(ani, block, sort_name, Show_fig)
                        rows = result[2].shape[0]
                        result_ft = pd.concat([df_header]*rows, ignore_index=True)
                        result_ft = pd.concat([result_ft, result[2]], axis = 1)
                        df_feature = pd.concat([df_feature, result_ft], ignore_index = True)
                    except:
                        print('No such path.')
                        
    return df_feature

                
            
            
        
        
    
    
    
    
    




    