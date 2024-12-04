import numpy as np
#=======================================
# Create PARSEC PARAMETER
#=======================================
R_LE       = 0.0155   #p1
X_UP       = 0.29663  #p2
Z_UP       = 0.06002  #p3
Z_XXUP     = -0.4515  #p4
X_LO       = 0.29663  #p5
Z_LO       = -0.06002 #p6
Z_XXLO     = 0.4514   #p7
Z_TE   	   = 0        #p8
delta_Z_TE = 0.0025   #P9
alfa_TE    = 0		  #p10
betha_TE   = 0.225    #p11
##Creates array of PARSEC parameters
PARSEC_var = np.array([R_LE,X_UP, Z_UP, Z_XXUP, X_LO,Z_LO,Z_XXLO, Z_TE, delta_Z_TE, alfa_TE, betha_TE])
#=======================================
# Create the number of x discrete points
#=======================================
x = np.linspace(0,1,100) #100 discrete points
#=======================================
# Create the X, C_up, C_low,b_up,b_low matrix
#=======================================
n=len(x)                    #number of discrete points
x_matrix = np.zeros ([n,6]) #creates nx6 matrix of zero
for i in  range (0,n):
    for j in range (0,6):
        x_matrix[i,j] = (x[i])**(j+0.5)  
######Create the C_up array
C_up = np.ones([6,6])       #creates 6x6 matrix contains number 1
for i in  range (1,6):
    C_up[5,i]=0
for i in  range (0,6):
    C_up[1,i] = (PARSEC_var[1]**(0.5+i))
    C_up[2,i] = 0.5+i
    C_up[3,i] = (PARSEC_var[1]**(-0.5+i))
    C_up[4,i] = (PARSEC_var[1]**(-1.5+i))
C_up[3,:] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]*C_up[3,:]
C_up[4,:] = [-1/4, 3/4, 15/4, 15/4, 63/4, 99/4]*C_up[4,:]
#####Create the C_low array
C_low = np.ones([6,6])
for i in  range (1,6):
    C_low[5,i]=0
for i in  range (0,6):
    C_low[1,i] = (PARSEC_var[4]**(0.5+i))
    C_low[2,i] = 0.5+i
    C_low[3,i] = (PARSEC_var[4]**(-0.5+i))
    C_low[4,i] = (PARSEC_var[4]**(-1.5+i))
C_low[3,:] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]*C_low[3,:]
C_low[4,:] = [-1/4, 3/4, 15/4, 15/4, 63/4, 99/4]*C_low[4,:]
#####Define b_up
b_up = np.zeros (6)
b_up[0] = PARSEC_var[7] +  PARSEC_var[8]/2
b_up[1] =  PARSEC_var[2]
b_up[2] = np.tan( PARSEC_var[9]- PARSEC_var[10]/2)
b_up[3] = 0
b_up[4] =  PARSEC_var[3]
b_up[5] = (2* PARSEC_var[0])**(0.5)
#####define b_low
b_low = np.zeros (6)
b_low[0] = PARSEC_var[7] -  PARSEC_var[8]/2
b_low[1] =  PARSEC_var[5]    
b_low[2] = np.tan( PARSEC_var[9] + PARSEC_var[10]/2)
b_low[3] = 0
b_low[4] =  PARSEC_var[6]   
b_low[5] = -(2* PARSEC_var[0])**(0.5)  
#calculate a_up and a_low
from numpy.linalg import inv,solve
a_up = solve(C_up,b_up)
a_low = solve(C_low,b_low)
#calculate Z_up and Z_low
z_up = np.dot(x_matrix,a_up)
z_low = np.dot(x_matrix,a_low)  
z_up = z_up.reshape([n,1])
z_low = z_low.reshape([n,1])
x = x.reshape([n,1])

#create a selig dat file format
print ("coordinate 1", '\n')
z_up = z_up.reshape([n,1])
z_low = z_low.reshape([n,1])
x = x.reshape([n,1])

x_reversed = np.flip (x)
z_up_reversed = np.flip(z_up)

coordinate_1 = np.append(x_reversed, z_up_reversed, axis=1) 
coordinate_2 = np.append(x, z_low, axis=1)
print (coordinate_1)

print ("coordinate 2",'\n')\

print (coordinate_2)

coordinate_final = np.append(coordinate_1[0:n-1,:], coordinate_2[0:n,:], axis=0)

print ('koordinat final','\n')
print (coordinate_final)

#export to selig txt file

np.savetxt('AIRFOIL_PARSEC.txt',coordinate_final, delimiter='  ')


#plot the data reference in a graph
import matplotlib.pyplot as plt
plt.plot(x, z_up, label='PARSEC Approaches')
plt.plot(x, z_low)
plt.grid(True)
plt.legend()
plt.xlabel('X ')
plt.ylabel('Z')
plt.show()

