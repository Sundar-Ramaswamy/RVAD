"""
            UTM Assignment
            Group number : 8
            Names (& CID) : Sundar Murugan Ramaswamy - murugan@net.chalmers.se
                            Gowtham Gunashekara - gowgun@net.chalmers.se
"""

#%% Imports

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%% Given
#General data
veh_speed= 80  #km/hr
altitude=1500 #m
rho_coolant= 1070 #kg/cu.m
cp_air= 1.004*1000 #J/KgK
rho_air= 1.205 #kg/cu.m
m_dot= (0.58*rho_coolant)/1000 #kg/s
Q_flow= 0.58*60 #L/min
v_flow=0.0834 #m/s
cp_coolant= 3.3*1000 #J/KgK
n=60 #number of cells
f=0.005
coolant_speed= 4000/60 #rps
L=4.3 #m
D=0.023 #m
frontal_area=0.489*0.98 #sq.m
k0_rad=4.58 #m^-4
k1_rad=40.11 #kgm-4s-1
k2_rad= -16.98 #kgm-1s-2
k0_cond=1.96 #m^-4
k1_cond=8.9 #kgm-4s-1
k2_cond= 6.37 #kgm-1s-2
R_air= 287 #gas constant
#Group specific data
Q_EM= 8000 #W
cell_heat_rate= 4 # W/cell
Q_dc= 6200 #W
Q_cond= 8000 #W
#%% Data maps
# Underhood fan




# Coolant pump
Dataset=pd.DataFrame(columns=['Q_flow','Pressure_Rise'])
Dataset.Q_flow=[0,25.2,34.8,44.12,65.32,93.83,120.67]
Dataset.Pressure_Rise=[127.92,130.47,np.nan,126.82,117.58,94.62,55.43]

#Task 1 - Temperature change across LT
#Pump
inlet_temp1=298 #degrees
outlet_temp1=298 #degrees %%since it is an adiabatic process

#battery
inlet_temp2= outlet_temp1 #deg
Q_b= cell_heat_rate*n #W (cumulative heat rate)
outlet_temp2= (Q_b/(m_dot*cp_coolant))+inlet_temp2
print("The outlet temperature at the battery is",outlet_temp2)

#Electric motor
inlet_temp3=outlet_temp2 #deg
outlet_temp3= (Q_EM/(m_dot*cp_coolant))+inlet_temp3
print("The outlet temperature at the EM is",outlet_temp3)

#DC-DC Converter
inlet_temp4=outlet_temp3 #deg
outlet_temp4= (Q_dc/(m_dot*cp_coolant))+inlet_temp4
print("The outlet temperature at the converter is",outlet_temp4)

Q_lt_rad= m_dot*cp_coolant*(inlet_temp1-outlet_temp4) #heat transferred across the LT Radiator
print("The heat transferred across the LT Radiator is",Q_lt_rad)

#Task 2 - Pressure drop across LT
#Pump
Dataset1=Dataset.interpolate(method='linear',limit_direction='forward') #code lines 51 to 53
print(Dataset1)
#battery
dP_batt= ((0.00375*np.power(Q_flow,2.5))-(0.2797*Q_flow))*1000
print("The pressure drop across the battery is",dP_batt)
#Electric motor
dP_EM= ((0.02*np.power(Q_flow,2))-(0.18*Q_flow))*1000
print("The pressure drop across the EM is",dP_EM)
#DC-DC Converter
dP_DCDC= ((0.001*np.power(Q_flow,3))-(0.0035*np.power(Q_flow,2))+(0.016*Q_flow))*1000
print("The pressure drop across the Converter is",dP_DCDC)
#frictional pressure drop
omega= 2*3.1416*coolant_speed
V= (D/2)*omega
delP_friction= (f*rho_coolant*np.power(v_flow,2)*L)/(2*D)
print("The frictional pressure drop is",delP_friction)
#whole radiator
dP_pump= 128.645*1000 #interpolated value pa
dP_rad=dP_pump-(dP_EM+dP_DCDC+dP_batt+delP_friction)
print("The pressure drop across the radiator is", dP_rad)


X = ['Pump','motor','DCDC','Battery','radiator','friction']
plt.figure()
plt.bar(range(len(X)),[dP_pump,dP_EM,dP_DCDC,dP_batt,dP_rad,delP_friction])
plt.xticks(range(len(X)),X)
plt.ylabel('change in pressure (Pa)')
plt.title('DelP in the system')

#Task 3 - air circuit
def air_circuit(m_air,Q_ht_rad):
    #before grille
    V= m_air/(0.42*rho_air*frontal_area)
    dyn_p_g= 0.5*rho_air*np.power(V,2)
    s_p_g=84300
    P_0= s_p_g+dyn_p_g
    T_0= 287.15 #k
    #after grille
    P_1= s_p_g-(0.35*dyn_p_g)  #due to loss in pressure when air passes through the grille which accounts to 35% of dynamic pressure
    T_1= T_0
    rho_1= P_1/(R_air*T_1)
    #after lt radiator
    dP_lt= (k0_rad*np.power(m_air,2))/rho_1 + (k1_rad*m_air)/rho_1 + k2_rad
    P_2= P_1-dP_lt
    T_2=(Q_lt_rad/(m_air*cp_air))+T_1
    rho_2=P_2/(R_air*T_2)
    #after condenser
    dP_cond= (k0_cond*np.power(m_air,2))/rho_2 + (k1_cond*m_air)/rho_2 + k2_cond
    P_3= P_2-dP_cond
    T_3=(Q_cond/(m_air*cp_air))+T_2
    rho_3=P_3/(R_air*T_3)
    #after ht radiator
    dP_ht= (k0_rad*np.power(m_air,2))/rho_3 + (k1_rad*m_air)/rho_3 + k2_rad
    P_4= P_3-dP_ht
    T_4=(Q_ht_rad/(m_air*cp_air))+T_3
    rho_4=P_4/(R_air*T_4)
    #after fan
    volFlowRate_4= m_air/rho_4
    #Dataset=pd.DataFrame(columns=['Q_flow','Pressure_Drop'])
    #Dataset.Q_flow=[0.0,1.6,2.0,3.0,4.0,4.9,5.05,5.2,6.0,7.2,8.0,9.9,11.1,12.3,13.2]
    #Dataset.Pressure_Drop=[4000.0,3500.0,3250.0,3000.0,2750.0,2500.0,2400.0,2300.0,2200.0,2100.0,2000.0,1500.0,1000.0,500.0,0.0]
    Q_f= np.array((0.0,1.6,2.0,3.0,4.0,4.9,5.05,5.2,6.0,7.2,8.0,9.9,11.1,12.3,13.2))
    P_D= np.array((4000.0,3500.0,3250.0,3000.0,2750.0,2500.0,2400.0,2300.0,2200.0,2100.0,2000.0,1500.0,1000.0,500.0,0.0))
    dP_fan= np.interp(volFlowRate_4,Q_f,P_D)
    P_5= P_4+dP_fan

    return dP_lt,dP_ht,dP_cond,P_1,P_2,P_3,P_4,P_5,T_2,T_3,T_4,rho_2,rho_3,rho_4,volFlowRate_4,dP_fan

#Task 4
#Design of experiment analysis
m_air = np.linspace(1.0, 12.0, num=50)
Q_ht_rad = np.linspace(40000.0 , 120000.0, num=50)
Data_set = pd.DataFrame(columns=['Q_ht_rad','mass_flowrate','dP_lt','dP_ht','dP_cond','P_1','P_2','P_3','P_4','P_5','T_2','T_3','T_4','rho_2','rho_3','rho_4','Q_4','dP_fan'])

#Loop over HT radiator heat rejected and mass flow rates
for x in Q_ht_rad:
    df = pd.DataFrame(columns=['mass_flowrate','dP_lt','dP_ht','dP_cond','P_1','P_2','P_3','P_4','P_5','T_2','T_3','T_4','rho_2','rho_3','rho_4','Q_4','dP_fan'])
    for y in m_air:
        #outputs = air_circuit()
        (dP_lt,dP_ht,dP_cond,P_1,P_2,P_3,P_4,P_5,T_2,T_3,T_4,rho_2,rho_3,rho_4,volFlowRate_4,dP_fan)=air_circuit(y,x)
        #Append datasets for loop 2
        df2 = pd.DataFrame([[x,y,dP_lt,dP_ht,dP_cond,P_1,P_2,P_3,P_4,P_5,T_2,T_3,T_4,rho_2,rho_3,rho_4,volFlowRate_4,dP_fan]]
                           ,columns=['Q_ht_rad','mass_flowrate','dP_lt','dP_ht','dP_cond','P_1','P_2','P_3','P_4','P_5','T_2','T_3','T_4','rho_2','rho_3','rho_4','Q_4','dP_fan'])
        df=df.append(df2)
    #append datasets for loop 1
    Data_set=Data_set.append(df)

#pressure contour plot
def air_circuit_1(m_air,Q_ht_rad):
    #before grille
    V= m_air/(0.42*rho_air*frontal_area)
    dyn_p_g= 0.5*rho_air*np.power(V,2)
    s_p_g=84300
    P_0= s_p_g+dyn_p_g
    T_0= 287.15 #k
    #after grille
    P_1= s_p_g-(0.35*dyn_p_g)  #due to loss in pressure when air passes through the grille which accounts to 35% of dynamic pressure
    T_1= T_0
    rho_1= P_1/(R_air*T_1)
    #after lt radiator
    dP_lt= (k0_rad*np.power(m_air,2))/rho_1 + (k1_rad*m_air)/rho_1 + k2_rad
    P_2= P_1-dP_lt
    T_2=(Q_lt_rad/(m_air*cp_air))+T_1
    rho_2=P_2/(R_air*T_2)
    #after condenser
    dP_cond= (k0_cond*np.power(m_air,2))/rho_2 + (k1_cond*m_air)/rho_2 + k2_cond
    P_3= P_2-dP_cond
    T_3=(Q_cond/(m_air*cp_air))+T_2
    rho_3=P_3/(R_air*T_3)
    #after ht radiator
    dP_ht= (k0_rad*np.power(m_air,2))/rho_3 + (k1_rad*m_air)/rho_3 + k2_rad
    P_4= P_3-dP_ht
    T_4=(Q_ht_rad/(m_air*cp_air))+T_3
    rho_4=P_4/(R_air*T_4)
    #after fan
    volFlowRate_4= m_air/rho_4
    #Dataset=pd.DataFrame(columns=['Q_flow','Pressure_Drop'])
    #Dataset.Q_flow=[0.0,1.6,2.0,3.0,4.0,4.9,5.05,5.2,6.0,7.2,8.0,9.9,11.1,12.3,13.2]
    #Dataset.Pressure_Drop=[4000.0,3500.0,3250.0,3000.0,2750.0,2500.0,2400.0,2300.0,2200.0,2100.0,2000.0,1500.0,1000.0,500.0,0.0]
    Q_f= np.array((0.0,1.6,2.0,3.0,4.0,4.9,5.05,5.2,6.0,7.2,8.0,9.9,11.1,12.3,13.2))
    P_D= np.array((4000.0,3500.0,3250.0,3000.0,2750.0,2500.0,2400.0,2300.0,2200.0,2100.0,2000.0,1500.0,1000.0,500.0,0.0))
    dP_fan= np.interp(volFlowRate_4,Q_f,P_D)
    P_5= P_4+dP_fan

    return P_5

Data_set1= pd.DataFrame(columns=['Q_ht_rad','mass_flowrate','P_5'])
A,B= np.meshgrid(m_air,Q_ht_rad)
lol= air_circuit_1(A,B)
plt.figure()
plt.contourf(A,B,lol,cmap='turbo');
plt.colorbar()
plt.xlabel('Mass flow rate (kg/s)')
plt.ylabel('Q_ht_rad (W)')
plt.title('Pressure downstream of the fan')

#temperature contour plot
def air_circuit_2(m_air,Q_ht_rad):
    #before grille
    V= m_air/(0.42*rho_air*frontal_area)
    dyn_p_g= 0.5*rho_air*np.power(V,2)
    s_p_g=84300
    P_0= s_p_g+dyn_p_g
    T_0= 287.15 #k
    #after grille
    P_1= s_p_g-(0.35*dyn_p_g)  #due to loss in pressure when air passes through the grille which accounts to 35% of dynamic pressure
    T_1= T_0
    rho_1= P_1/(R_air*T_1)
    #after lt radiator
    dP_lt= (k0_rad*np.power(m_air,2))/rho_1 + (k1_rad*m_air)/rho_1 + k2_rad
    P_2= P_1-dP_lt
    T_2=(Q_lt_rad/(m_air*cp_air))+T_1
    rho_2=P_2/(R_air*T_2)
    #after condenser
    dP_cond= (k0_cond*np.power(m_air,2))/rho_2 + (k1_cond*m_air)/rho_2 + k2_cond
    P_3= P_2-dP_cond
    T_3=(Q_cond/(m_air*cp_air))+T_2
    rho_3=P_3/(R_air*T_3)
    #after ht radiator
    dP_ht= (k0_rad*np.power(m_air,2))/rho_3 + (k1_rad*m_air)/rho_3 + k2_rad
    P_4= P_3-dP_ht
    T_4=(Q_ht_rad/(m_air*cp_air))+T_3
    rho_4=P_4/(R_air*T_4)
    #after fan
    volFlowRate_4= m_air/rho_4
    #Dataset=pd.DataFrame(columns=['Q_flow','Pressure_Drop'])
    #Dataset.Q_flow=[0.0,1.6,2.0,3.0,4.0,4.9,5.05,5.2,6.0,7.2,8.0,9.9,11.1,12.3,13.2]
    #Dataset.Pressure_Drop=[4000.0,3500.0,3250.0,3000.0,2750.0,2500.0,2400.0,2300.0,2200.0,2100.0,2000.0,1500.0,1000.0,500.0,0.0]
    Q_f= np.array((0.0,1.6,2.0,3.0,4.0,4.9,5.05,5.2,6.0,7.2,8.0,9.9,11.1,12.3,13.2))
    P_D= np.array((4000.0,3500.0,3250.0,3000.0,2750.0,2500.0,2400.0,2300.0,2200.0,2100.0,2000.0,1500.0,1000.0,500.0,0.0))
    dP_fan= np.interp(volFlowRate_4,Q_f,P_D)
    P_5= P_4+dP_fan

    return T_4

Data_set1= pd.DataFrame(columns=['Q_ht_rad','mass_flowrate','T_5'])
A,B= np.meshgrid(m_air,Q_ht_rad)
lol= air_circuit_2(A,B)
plt.figure()
plt.contourf(A,B,lol,cmap='turbo');
plt.colorbar()
plt.xlabel('Mass flow rate (kg/s)')
plt.ylabel('Q_ht_rad (W)')
plt.title('temperature downstream of the fan')

#density contour plot
def air_circuit_3(m_air,Q_ht_rad):
    #before grille
    V= m_air/(0.42*rho_air*frontal_area)
    dyn_p_g= 0.5*rho_air*np.power(V,2)
    s_p_g=84300
    P_0= s_p_g+dyn_p_g
    T_0= 287.15 #k
    #after grille
    P_1= s_p_g-(0.35*dyn_p_g)  #due to loss in pressure when air passes through the grille which accounts to 35% of dynamic pressure
    T_1= T_0
    rho_1= P_1/(R_air*T_1)
    #after lt radiator
    dP_lt= (k0_rad*np.power(m_air,2))/rho_1 + (k1_rad*m_air)/rho_1 + k2_rad
    P_2= P_1-dP_lt
    T_2=(Q_lt_rad/(m_air*cp_air))+T_1
    rho_2=P_2/(R_air*T_2)
    #after condenser
    dP_cond= (k0_cond*np.power(m_air,2))/rho_2 + (k1_cond*m_air)/rho_2 + k2_cond
    P_3= P_2-dP_cond
    T_3=(Q_cond/(m_air*cp_air))+T_2
    rho_3=P_3/(R_air*T_3)
    #after ht radiator
    dP_ht= (k0_rad*np.power(m_air,2))/rho_3 + (k1_rad*m_air)/rho_3 + k2_rad
    P_4= P_3-dP_ht
    T_4=(Q_ht_rad/(m_air*cp_air))+T_3
    rho_4=P_4/(R_air*T_4)
    #after fan
    volFlowRate_4= m_air/rho_4
    #Dataset=pd.DataFrame(columns=['Q_flow','Pressure_Drop'])
    #Dataset.Q_flow=[0.0,1.6,2.0,3.0,4.0,4.9,5.05,5.2,6.0,7.2,8.0,9.9,11.1,12.3,13.2]
    #Dataset.Pressure_Drop=[4000.0,3500.0,3250.0,3000.0,2750.0,2500.0,2400.0,2300.0,2200.0,2100.0,2000.0,1500.0,1000.0,500.0,0.0]
    Q_f= np.array((0.0,1.6,2.0,3.0,4.0,4.9,5.05,5.2,6.0,7.2,8.0,9.9,11.1,12.3,13.2))
    P_D= np.array((4000.0,3500.0,3250.0,3000.0,2750.0,2500.0,2400.0,2300.0,2200.0,2100.0,2000.0,1500.0,1000.0,500.0,0.0))
    dP_fan= np.interp(volFlowRate_4,Q_f,P_D)
    P_5= P_4+dP_fan

    return rho_4

Data_set1= pd.DataFrame(columns=['Q_ht_rad','mass_flowrate','rho_4'])
A,B= np.meshgrid(m_air,Q_ht_rad)
lol= air_circuit_3(A,B)
plt.figure()
plt.contourf(A,B,lol,cmap='turbo');
plt.colorbar()
plt.xlabel('Mass flow rate (kg/s)')
plt.ylabel('Q_ht_rad (W)')
plt.title('density downstream of the fan')

#task 5
m_air = np.linspace(1.0, 12.0,120)
Q_ht_rad = np.linspace(40000,120000,120)
Data_set4 = pd.DataFrame(columns=['Q_ht_rad','mass_flowrate','dP_lt','dP_ht','dP_cond','P_1','P_2','P_3','P_4','P_5','T_2','T_3','T_4','rho_2','rho_3','rho_4','Q_4','dP_fan'])

#Loop over HT radiator heat rejected and mass flow rates
for x in Q_ht_rad:
    df = pd.DataFrame(columns=['mass_flowrate','dP_lt','dP_ht','dP_cond','P_1','P_2','P_3','P_4','P_5','T_2','T_3','T_4','rho_2','rho_3','rho_4','Q_4','dP_fan'])
    for y in m_air:
        #outputs = air_circuit()
        (dP_lt,dP_ht,dP_cond,P_1,P_2,P_3,P_4,P_5,T_2,T_3,T_4,rho_2,rho_3,rho_4,volFlowRate_4,dP_fan)=air_circuit(y,x)
        #Append datasets for loop 2
        df2 = pd.DataFrame([[x,y,dP_lt,dP_ht,dP_cond,P_1,P_2,P_3,P_4,P_5,T_2,T_3,T_4,rho_2,rho_3,rho_4,volFlowRate_4,dP_fan]]
                           ,columns=['Q_ht_rad','mass_flowrate','dP_lt','dP_ht','dP_cond','P_1','P_2','P_3','P_4','P_5','T_2','T_3','T_4','rho_2','rho_3','rho_4','Q_4','dP_fan'])
        df=df.append(df2)
    #append datasets for loop 1
    Data_set4=Data_set4.append(df)

P_0=84300
x=0
m_1=np.zeros(120)

for a in Q_ht_rad:
    y= Data_set4[Data_set4.Q_ht_rad==a].P_5.values
    T= Data_set4[Data_set4.Q_ht_rad==a].mass_flowrate.values
    for b in range(0,120):
        if y[b] < P_0:
            y1 = y[b-1]
            y2 = y[b]
            m1 = T[b-1]
            m2 = T[b]

            m_1[x] = ((P_0 - y2) * ((m2 - m1) / (y2 - y1))) + m2
            x=x+1
            break

plt.figure()
plt.grid()
plt.plot(m_1,Q_ht_rad,'x-')
plt.ylabel('Q_ht_rad (W)')
plt.xlabel('Mass flow rate (kg/s)')
plt.title('Mass Flow Rate vs Heat Transferred')
