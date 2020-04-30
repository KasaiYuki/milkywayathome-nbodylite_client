import numpy as np
import math as m
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

Z_0 = 27 #pc
h_d = 250 #pc
l_d = 2250 #pc
l_td = 3500 #pc
h_td = 700 #pc
p_0td = 8 #%
p_0halo = .125 #%
q = .8
n = 3.5

A_n = 70 #pc
A_s = 170 #pc
lambda_n = 6.3 #kpc
lambda_s = 7.8 #kpc
r_s = 14 #kpc
a_0 = 1000 #pc
p_0 = .14
def p_d(r,z):
    return np.exp(-r/l_d) * np.exp(-abs(z-z_w(r))/h_d)
def p_td(r,z):
    return np.exp(-r/l_td)* np.exp(-abs(z-z_w(r))/h_td)
def p_s(r,z):
    d = (r**2 + (z/.8)**2)**.5
    return 1/(1000**3.5 + d**3.5)
def z_w(r):
    r_kpc = r / 1000.
    if r_kpc < 12.05:
        return 70*np.sin(2*np.pi*(r_kpc-2.225)/6.3)
    else:
        return 170*np.sin(2*np.pi*(r_kpc-.35)/7.8)

def find_density(r,z):
    p_0 = .14
    z = 0
    return p_0*(0.91875 * p_d(r,z) + .08 * p_td(r,z) + .00125* p_s(r,z))

def dist(mag): #returns in PC
    exponent = (mag - 4.2)/5.
    lhs = 10**exponent
    lhs = float(lhs) * 10.
    return lhs


#phi is the angle from the +z axis going down [0, pi]
#theta is the angle from the +x axis going around [0, 2*pi]
#p is radius (r) [0, inf)

def numerical_int(r_start, r_end, b_start, b_end, l_start, l_end):
    phi_start = np.radians(90. - b_start)
    phi_end = np.radians(90. - b_end)
    theta_start = np.radians(l_start)
    theta_end = np.radians(l_end)
    dr = (r_end-r_start)/10.
    dphi = (phi_end-phi_start)/10.
    dtheta =(theta_end-theta_start)/10.
    radius = np.arange(r_start, r_end, dr)
    phi = np.arange(phi_start, phi_end, dphi)
    theta = np.arange(theta_start, theta_end, dtheta)
        
    volume = 0
    total = 0.
    for i in radius:
        for j in theta:
            for k in phi:
                x_sun = i*np.cos(j)*np.sin(k)
                y_sun = i*np.sin(j)*np.sin(k)
                z_sun = i*   1     *np.cos(k)          
                #shifted_R = ((x_sun-8000)**2 + (y_sun)**2)**.5
                shifted_R = ((x_sun-0)**2 + (y_sun)**2)**.5
                dV = i**2 * np.sin(k) *dr*dphi*dtheta
                volume += dV
                total  += find_density(shifted_R,z_sun) * dV 
    if is_starcounts == "1":
        return total/1
    else:
        return total/volume

if __name__ == "__main__":
    print("Working on Disk Model.\nWant linear scale (0) or log scale (1)?")
    is_log = "1"
    #is_log = input()
    print("Want density (0) or star counts (1)?")
    is_starcounts = "0"
    #is_starcounts = input() 
    magnitude = [i/2 +16 for i in range(0,15)]
    l = [i*5 for i in range(0, 72)]
    b = [i*5 +2.5 for i in range(1,6)] + [i*-5 -2.5 for i in range(1,6)]
    fp = sys.stdout
    fp.flush()
    t1 = time.time()
    for inclination in b[0:1]:
        b_low = inclination - 2.5
        b_high = inclination + 2.5
        data = []
        temp = []
        for m in magnitude:
            r_low = dist(m - .25)
            r_high = dist(m + .25)
            
            for angle in l:
                l_low = angle
                l_high = angle + 5
                
                density = numerical_int(r_low, r_high, b_high, b_low, l_low, l_high)
                data.append((m, angle, density))
            temp.append((m, angle, density))
        
        r = [radius for (radius, a, b) in data]
        theta = [np.radians(angle) for (a, angle, b) in data]
        colors =[shade for (a, b, shade) in data]    
        
        area = 10
        fig = plt.figure()
        
        ax = fig.add_subplot(111, projection='polar')
        ax.set_ylim(15,23)
        index = 0
        
        #print(colors, "\n")
        
        if (is_log == "1"):
            c = ax.scatter(theta, r, c=colors, s=area, cmap='inferno_r', alpha=1, norm = matplotlib.colors.LogNorm())
        else:        
            c = ax.scatter(theta, r, c=colors, s=area, cmap='inferno_r', alpha=1)
        plt.colorbar(c)
        print("Predicted Stellar Density for b " + str(inclination) + ".png")
        
        fig_title = "Predicted Disk Model for b " + str(inclination)
        if is_starcounts != "0":
            fig_title +=  "  Star Counts "
        else:
            fig_title += " Density "
        if is_log != "0":
            fig_title += " Logscale"
        else:
            fig_title += " LinearScale"  
        fig.suptitle(fig_title)
        fig.savefig(fig_title+  ".png")
        plt.close()    
    t2 = time.time()
    print(t2-t1)

    