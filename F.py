#reading the dataframe
import pandas as pd
print('Enter the drive')
drive = input()
print('Enter the folder name:')
folder= input()
str_folder = str(folder)+'_without_result'+'.csv'
df_save=pd.read_csv(str_folder)
df_save['scf_cc']=0.000
df_save['scf_bc']=0.000
df_save['scf_bs']=0.000
df_save['scf_cs']=0.000
import os
import math
#result file scrapping
for each in range(len(df_save)):
    alpha = df_save['alpha'][each]
    beta =df_save['beta'][each]
    gamma = df_save['gamma'][each]
    tau = df_save['tau'][each]
    C_D = df_save['C_D'][each]
    #C_L = 0.5*alpha*C_D
    B_D = C_D*beta
    C_T = 0.5*C_D/gamma
    B_T = tau*C_T 
    C_R = C_D/2
    B_R = B_D/2
    ao = 0.2*math.sqrt(B_R*B_T)
    lrmin = max(4,0.4*C_T )
    if ao<=lrmin:
      a=lrmin
    else:
      a=ao
    lrmin = max(4,0.4*B_T )
    if ao<=lrmin:
      a_b=lrmin
    else:
      a_b=ao
    b_b = 0.65*math.sqrt(B_R*B_T)
    b_cc = 0.4*math.sqrt(math.sqrt(B_R*C_R*B_T*C_T))
    b_cs = min(0.09*C_R,math.pi*C_R/36)
    if b_b <= a_b + 0.6*B_T:
      b_b = a_b + 0.6*B_T
    if b_cc <= a + 0.6*B_T:
      b_cc = a + 0.6*B_T  
    if b_cs <= a + 0.6*B_T:
      b_cs = a+ 0.6*B_T   
    a_b =0.4*C_T 
    a = 0.4*C_T
    b_b = 1.4*C_T
    b_cc = b_b
    b_cs = b_b  
    inp_dir = each
    path = ''''''+str(drive.upper()) +''':\\'''+str(folder)+'\\'+str(inp_dir)+'/'
    directory = ''''''+str(drive.upper()) +''':\\'''+str(folder)+'\\'
    file1 = open(path+'nodeindex.txt', 'r')
    Lines = file1.readlines()
    df = pd.DataFrame(columns=['node','x','y','z'])
    for i in Lines:
        out =i.split()
        if "-" in out[0]:
            j=out[0].split("-")
            out.pop(0)
            k =j[1]+'-'+j[2]
            j.pop(1)
            j.pop(1)
            j.insert(1,k)
            out.insert(0,j[0])
            out.insert(1,j[1])
        if len(out)==4:
            df = df.append(pd.Series(out, index=df.columns), ignore_index=True)
        elif len(out)==3:
            out.insert(3,0)
            df = df.append(pd.Series(out, index=df.columns), ignore_index=True)
    #converting object type into float
    df['x'] = df['x'].astype(float)
    df['y'] = df['y'].astype(float)
    df['z'] = df['z'].astype(float)
    #Scrapping
    df.loc[df['x']<10**-5, 'x']=0
    #Indexing 
    list_z,list_y,list_zx,list_zy = [],[],[],[]
    index =[]
    for i in range(len(df['x'])):
        if df['x'][i] == 0:
            list_z.append(df['z'][i])
            list_y.append(df['y'][i]) 
        else:
            if df['z'][i] == 0:
                list_zx.append(df['x'][i])
                list_zy.append(df['y'][i])
    list_z.sort()
    list_y.sort()
    list_zx.sort()
    list_zy.sort()
    for i in range(len(df['z'])):
        if df['z'][i] == list_z[-1]: index.append('a')
        elif df['z'][i] == list_z[-2]: index.append('b')
        elif df['z'][i] == list_z[-3]: index.append('c')
        else: 
            if df['y'][i]==list_y[-1]: index.append('d')
            elif df['y'][i]==list_y[-2]: index.append('e')
            elif df['y'][i]==list_y[-3]: index.append('f')

            else:
                if df['y'][i] == list_zy[-1]: index.append('g')
                elif df['y'][i] == list_zy[-2]: index.append('h')
                elif df['y'][i] == list_zy[-3]: index.append('i')
                else: 
                    if df['x'][i]==list_zx[-1]: index.append('j')
                    elif df['x'][i]==list_zx[-2]: index.append('k')
                    elif df['x'][i]==list_zx[-3]: index.append('l')
    df['index'] = index
    sorted_df =df.sort_values('index').reset_index(drop=True)

    file2 = open(path +'result.txt', 'r')   
    lines = file2.readlines()
    stress_x=[]
    stress_y=[]
    stress_z=[]
    stress_xy=[]
    count =0
    for line in lines[1:19]:
        p=line.split(",")
        p.append(p[2].strip('\n' or '\t'))
        p.pop(2)
        if count == 0 or count== 1 or count== 2:
            stress_z.append(p[-1])
            stress_y.append(0)
            stress_x.append(0)
            stress_xy.append(0)
        elif count == 3 or count== 4 or count== 5:
            stress_y.append(p[-1])
            stress_x.append(0)
            stress_z.append(0)
            stress_xy.append(0)
        elif count == 6 or count== 7 or count== 8:
            stress_y.append(p[-1])
            stress_x.append(0)
            stress_z.append(0)
            stress_xy.append(0)
        elif count==9 or count== 12 or count== 15:
            stress_x.append(p[-1])
            stress_z.append(0)
        elif count==10 or count== 13 or count== 16:
            stress_xy.append(p[-1])
        else:
            stress_y.append(p[-1])
        count+=1
    #clean data in x,y,z
    stress_x=[float(i) for i in stress_x]
    stress_y=[float(i) for i in stress_y]
    stress_z=[float(i) for i in stress_z]
    stress_xy=[float(i) for i in stress_xy]
    sorted_df['stress_x'] = stress_x
    sorted_df['stress_y'] = stress_y
    sorted_df['stress_z'] = stress_z
    sorted_df['stress_xy'] = stress_xy
    #stress extrapolation
    #chord_crown
    temp_b = (4*(sorted_df['stress_z'][1])-sorted_df['stress_z'][2]-(3*sorted_df['stress_z'][0]))/(2*(b_cc-a))
    temp_a = ((sorted_df['stress_z'][1]-sorted_df['stress_z'][0])-(temp_b*(b_cc-a)))/((b_cc-a)**2)
    SCF_chord_crown = abs((temp_a*(((2*b_cc)-a)**2))+(temp_b*((2*b_cc)-a))+sorted_df['stress_z'][0])
    #brace_crown
    temp_b = (4*(sorted_df['stress_y'][4])-sorted_df['stress_y'][5]-(3*sorted_df['stress_y'][3]))/(2*(b_b-a_b))
    temp_a = ((sorted_df['stress_y'][4]-sorted_df['stress_y'][3])-(temp_b*(b_b-a_b)))/((b_b-a_b)**2)
    SCF_brace_crown = abs((temp_a*(((2*b_b)-a_b)**2))+(temp_b*((2*b_b)-a_b))+sorted_df['stress_y'][3])
    #brace_saddle
    temp_b = (4*(sorted_df['stress_y'][7])-sorted_df['stress_y'][8]-(3*sorted_df['stress_y'][6]))/(2*(b_b-a_b))
    temp_a = ((sorted_df['stress_y'][7]-sorted_df['stress_y'][6])-(temp_b*(b_b-a_b)))/((b_b-a_b)**2)
    SCF_brace_saddle = abs((temp_a*(((2*b_b)-a_b)**2))+(temp_b*((2*b_b)-a_b))+sorted_df['stress_y'][6])
    #chord_saddle
    stress_theta =[]
    for i in range(9,12,1):
        theta=math.atan(sorted_df['y'][i]/sorted_df['x'][i])
        stress_theta.append(((sorted_df['stress_x'][i]+sorted_df['stress_y'][i])*0.5)-(((sorted_df['stress_x'][i]-sorted_df['stress_y'][i])*0.5)*math.cos(2*theta))-(sorted_df['stress_xy'][i]*math.sin(2*theta)))
    temp_b = (4*(stress_theta[1])-stress_theta[2]-(3*stress_theta[0]))/(2*(b_cs-a))
    temp_a = ((stress_theta[1]-stress_theta[0])-(temp_b*(b_cs-a)))/((b_cs-a)**2)
    SCF_chord_saddle = abs((temp_a*(((2*b_cs)-a)**2))+(temp_b*((2*b_cs)-a))+stress_theta[0])
    SCF_quadratic = [SCF_chord_crown,SCF_brace_crown,SCF_brace_saddle,SCF_chord_saddle]
    #linear SCF
    SCF_chord_crown = abs(sorted_df['stress_z'][2]+(a*(sorted_df['stress_z'][2]-sorted_df['stress_z'][1])/(b_cc-a)))
    SCF_brace_crown = abs(sorted_df['stress_y'][5]+(a_b*(sorted_df['stress_y'][5]-sorted_df['stress_y'][4])/(b_b-a)))
    SCF_brace_saddle = abs(sorted_df['stress_y'][8]+(a_b*(sorted_df['stress_y'][8]-sorted_df['stress_y'][7])/(b_b-a)))
    SCF_chord_saddle = abs(stress_theta[2]+(a*(stress_theta[2]-stress_theta[1])/(b_cs-a)))
    SCF_linear = [SCF_chord_crown,SCF_brace_crown,SCF_brace_saddle,SCF_chord_saddle]
    #Appending SCF
    SCF =[]
    import numpy as np
    if np.sign(sorted_df['stress_z'][0])==np.sign(sorted_df['stress_z'][1]) and np.sign(sorted_df['stress_z'][0])==np.sign(sorted_df['stress_z'][2]):
        SCF.append(SCF_quadratic[0])
    else:
        SCF.append(SCF_quadratic[0])
    if np.sign(sorted_df['stress_y'][3])==np.sign(sorted_df['stress_y'][4]) and np.sign(sorted_df['stress_y'][3])==np.sign(sorted_df['stress_y'][5]):
        SCF.append(SCF_quadratic[1])
    else:
        SCF.append(SCF_quadratic[1])  
    if np.sign(sorted_df['stress_y'][6])==np.sign(sorted_df['stress_y'][7]) and np.sign(sorted_df['stress_y'][6])==np.sign(sorted_df['stress_y'][8]):
        SCF.append(SCF_quadratic[2])
    else:
        SCF.append(SCF_quadratic[2])
    if np.sign(stress_theta[0])==np.sign(stress_theta[1]) and np.sign(stress_theta[0])==np.sign(stress_theta[2]):
        SCF.append(SCF_quadratic[3])
    else:
        SCF.append(SCF_quadratic[3])
    # adding results to the df 
    df_save['scf_cc'][each]=SCF[0]
    df_save['scf_bc'][each]=SCF[1]
    df_save['scf_bs'][each]=SCF[2]
    df_save['scf_cs'][each]=SCF[3]
# saving the dataframe
str_folder = str(folder)+'_with_result'+'.csv'
df_save.to_csv(str_folder)
