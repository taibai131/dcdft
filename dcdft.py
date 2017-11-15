import numpy as np
import matplotlib.pyplot as plt 

# Open data
t, dB = np.loadtxt("dtvir.data", unpack=True)
avgdB = np.mean(dB)
f = dB - avgdB		# selisih magnitudo dengan magnitudo rata-rata

# Jumlah data
n = len(t)

# Plot dB vs JD
plt.figure(0)
plt.plot(t, dB, 'k.')
plt.title("$\Delta$B vs JD 2441400+")

# Find frequency Nyquist
diff_t = np.diff(t)	# Selisih setiap waktu pengamatan
wa = 1/(np.min(diff_t)*2) # w = 1/2*Pmin

# Selang 
s = .0001 			# Selang/Step, semakin kecil semakin baik, namun tidak boleh terlalu kecil
dt = t[n-1]-t[0] 	# Waktu pengamatan, dalam satuan hari

# metode Date Compensated Discrete Fourier Transform (DCDFT)
a0 = np.sqrt(1/n)
S = []
H = []
C = []
alpha = ((2.*(n-3)*dt*wa)/(3.*(n-4)))

ws = np.arange(s,wa,s)
for w in ws:
	sum_cos2x = sum_cosx = sum_cosxsinx = sum_sin2x = sum_sinx = sum_f_cos = sum_f_sin = sum_f2 = 0
	for i in range(n):
		x = 2.*np.pi*w*t[i]
		sum_cos2x += np.cos(x)**2
		sum_cosx += np.cos(x)
		sum_cosxsinx += (np.cos(x)*np.sin(x))
		sum_sin2x += np.sin(x)**2
		sum_sinx += np.sin(x)
		sum_f_cos += f[i]*np.cos(x)
		sum_f_sin += f[i]*np.sin(x)
		sum_f2 += f[i]**2
		pass
	M = sum_cosxsinx - (a0**2)*(sum_sinx)*(sum_cosx)
	a1 = np.sqrt(1./(sum_cos2x - (a0**2)*(sum_cosx**2)))
	a2 = np.sqrt(1./(sum_sin2x - (a0**2)*(sum_sinx**2) - (a1**2)*M))
	c1 = a1*sum_f_cos
	c2 = a2*sum_f_sin - (a1*a2*c1*M)
	S_ = ( c1**2. + c2**2. )/sum_f2
	G = (-1.*((n-3)/2.)*np.log(1.-(S_)))
	H_ = ((n-4)*1./(n-3))*(G + np.exp(-1.*G)-1.)
	S.append(S_)
	H.append(H_)
	C.append(100*(1-np.exp(-H_))**alpha)

fig, ax1 = plt.subplots()
ax1.plot(np.arange(s,wa,s),H,"b-", label="H")
ax1.plot(np.arange(s,wa,s),S,"r-", label="Spectral Correlation")
ax2 = ax1.twinx()
ax2.plot(np.arange(s,wa,s),C,"k-.", lw=1.2, label="Confident Level (%)")
fig.tight_layout()
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")
ax1.set_xlabel('$frequency$')
ax1.set_ylabel('H($\omega$)')
ax1.set_xlim(0,wa)
ax2.set_ylabel('Confident Level (%)')
ax2.grid(color='k', linestyle='--', linewidth=.5)
ax1.set_title("Modified Periodogram or graph of the function H($\omega$). The confidence levels are also shown.")

# Mencari Periode dari 
P = 1/ws[np.argmax(H)]

fase_ = (t-t[0])/P
fase = (fase_) - np.floor(fase_)
fase = fase - fase[np.argmax(dB)]
fase_negatif = fase[fase < 0.]
fase = fase.tolist()
arg_fase_negatif = []
for neg in fase_negatif:
	# arg_fase_negatif.append(fase.index(neg))
	index = fase.index(neg)
	fase[index] = 1 + neg

plt.figure(2)
plt.plot(fase,dB,'b.')
plt.xlabel('fase') 
plt.xlim(0,1)
plt.ylabel('$\Delta$B')
plt.title('Kurva Cahaya DT Vir')

## Printing hasil perhitungan
# Print Periode
print ("Frequency Nyquist = {:.3f}".format(wa))
print ("Periode dari bintang DT Vir adalah {:.3f}".format(P))
# Print data
print ("\nTabel 1. Data DT Vir")
print ("============================================")
print ("t\t\tdelta_B\t\tf\tfase")
for i in range(n):
	print ("{}\t{}\t{:.3f}\t{:.3f}".format(t[i],dB[i],f[i],fase[i]))

# Show plot
plt.show()