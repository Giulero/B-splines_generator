from bspline import *

ctr = np.array([(1.925000,-1.250000),
(1.800000,-1.175000),
(1.700000,-1.099999),
(1.575000,-1.025000),
(1.475000,-0.950000),
(1.350000,-0.875000),
(1.250000,-0.800000),
(1.125000,-0.724999),
(1.025001,-0.650000),
(1.000000,-0.500000),
(0.925000,-0.375000),
(0.875000,-0.250000),
(0.825000,-0.125000),
(0.775001,0.000000),
(0.800000,0.125000),
(0.800000,0.250000),
(0.800000,0.375000),
(0.825000,0.500000),
(0.825000,0.625000),
(0.825000,0.750000),
(0.825000,0.875000),
(0.850000,1.000000),
(0.850000,1.075000),
(0.775001,1.200000),
(0.725000,1.325000),
(0.675000,1.450000),
(0.625000,1.575000),
(0.600000,1.700000),
(0.600000,1.825000),
(0.600000,1.950000),
(0.575000,2.075000),
(0.575000,2.200000),
(0.575000,2.225000),
(0.625000,2.350000),
(0.650001,2.475000),
(0.700000,2.600000),
(0.750000,2.725000),
(0.800000,2.850000),
(0.850000,2.975000),
(0.925000,3.075000),
(1.025001,3.000000),
(1.150001,2.925000),
(1.275001,2.850000),
(1.350000,2.800000)]).transpose()

x=ctr[0,:]
y=ctr[1,:]

m = len(x)
p = 3; #polygrade

s = np.linspace(0, m, m+p+1-2*p)

# Knot vector
knots = np.zeros(m+p+1)
knots[:(p+1)]=s[0]
knots[-(p+1):]=s[-1]
knots[(p+1):-(p+1)]= s[1:-1]

u3=np.linspace(0,m,(max(m*2,70)))

# Vectors used to evaluate the b-splines
out = np.zeros((2, len(u3)))
out_der = np.zeros((2, len(u3)))

for i in range(len(u3)):
	for j in (range(m)):
		out[0,i] += ctr[0,j]*BasisFunc(u3[i], j, knots, p)
		out[1,i] += ctr[1,j]*BasisFunc(u3[i], j, knots, p)

for i in range(len(u3)):
	for j in (range(m)):
		out_der[0,i] += ctr[0,j]*BasisFunc_der(u3[i], j, knots, p)
		out_der[1,i] += ctr[1,j]*BasisFunc_der(u3[i], j, knots, p)

print('out_der', out_der)

plt.plot(out[0],out[1],'b',linewidth=2.0,label='B-spline curve')
plt.plot(x,y,'k--',label='Control polygon',marker='o',markerfacecolor='red')

plt.savefig('bspline.pdf')
plt.show()
plt.plot(out_der[0],'b',linewidth=2.0,label='B-spline curve')
plt.savefig('v_x.pdf')
plt.show()
plt.plot(out_der[1],'b',linewidth=2.0,label='B-spline curve')
plt.savefig('v_y.pdf')
plt.show()
