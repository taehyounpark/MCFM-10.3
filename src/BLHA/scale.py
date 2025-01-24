
import math as m

class Vec4:

    def __init__(self,E=0,px=0,py=0,pz=0):
        self.E = E
        self.px = px
        self.py = py
        self.pz = pz

    def __getitem__(self,i):
        if i == 0: return self.E
        if i == 1: return self.px
        if i == 2: return self.py
        if i == 3: return self.pz
        raise Exception('Vec4D')

    def __setitem__(self,i,v):
        if i == 0: self.E = v
        if i == 1: self.px = v
        if i == 2: self.py = v
        if i == 3: self.pz = v

    def __repr__(self):
        return '({0},{1},{2},{3})'.format(self.E,self.px,self.py,self.pz)

    def __str__(self):
        return '({0},{1},{2},{3})'.format(self.E,self.px,self.py,self.pz)

    def __add__(self,v):
        return Vec4(self.E+v.E,self.px+v.px,self.py+v.py,self.pz+v.pz)

    def __sub__(self,v):
        return Vec4(self.E-v.E,self.px-v.px,self.py-v.py,self.pz-v.pz)

    def __neg__(self):
        return Vec4(-self.E,-self.px,-self.py,-self.pz)

    def __mul__(self,v):
        if isinstance(v,Vec4):
            return self.E*v.E-self.px*v.px-self.py*v.py-self.pz*v.pz
        return Vec4(self.E*v,self.px*v,self.py*v,self.pz*v)

    def __rmul__(self,v):
        if isinstance(v,Vec4):
            return self.E*v.E-self.px*v.px-self.py*v.py-self.pz*v.pz
        return Vec4(self.E*v,self.px*v,self.py*v,self.pz*v)

    def __div__(self,v):
        return Vec4(self.E/v,self.px/v,self.py/v,self.pz/v)

    def Vec3(self):
        return Vec4(0.0,self.px,self.py,self.pz)

    def M2(self):
        return self*self

    def M(self):
        m2 = self.M2()
        return m.sqrt(m2 if m2>0.0 else -m2)

    def P2(self):
        return self.px*self.px+self.py*self.py+self.pz*self.pz

    def P(self):
        return m.sqrt(self.P2())

    def PT2(self):
        return self.px*self.px+self.py*self.py

    def PT(self):
        return m.sqrt(self.PT2())

    def Theta(self) :
        return m.acos(self.pz/self.P())

    def Phi(self) :
        if self.px==0 and self.py==0:
            return 0.0
        else:
            return m.atan2(self.py,self.px)

    def Cross(self,v):
        return Vec4(0.0,
                    self.py*v.pz-self.pz*v.py,
                    self.pz*v.px-self.px*v.pz,
                    self.px*v.py-self.py*v.px)

    def Boost(self,v):
        rsq = self.M()
        v0 = (self.E*v.E-self.px*v.px-self.py*v.py-self.pz*v.pz)/rsq;
        c1 = (v.E+v0)/(rsq+self.E)
        return Vec4(v0,
                    v.px-c1*self.px,
                    v.py-c1*self.py,
                    v.pz-c1*self.pz)

    def BoostBack(self,v):
        rsq = self.M()
        v0 = (self.E*v.E+self.px*v.px+self.py*v.py+self.pz*v.pz)/rsq;
        c1 = (v.E+v0)/(rsq+self.E)
        return Vec4(v0,
                    v.px+c1*self.px,
                    v.py+c1*self.py,
                    v.pz+c1*self.pz)

class Rotation:
    def __init__(self,v1,v2):
        a = v1.Vec3()/v1.P()
        b = v2.Vec3()/v2.P()
        self.x = a
        self.y = b+a*(a*b)
        if self.y.P2() != 0.:
            self.y = self.y/self.y.P()
        self.ct = -self.x*b
        self.st = -self.y*b
    def __mul__(self,v):
        vx = -self.x*v
        vy = -self.y*v
        vt = v-self.x*vx-self.y*vy
        return vt+(self.ct*vx-self.st*vy)*self.x \
            +(self.st*vx+self.ct*vy)*self.y

def sign(x):
    if x<0.0: return -1
    return 1

def LT(a,b,c):
    t = a[1]*b[2]*c[3]+a[2]*b[3]*c[1]+a[3]*b[1]*c[2] \
       -a[1]*b[3]*c[2]-a[3]*b[2]*c[1]-a[2]*b[1]*c[3]
    x = -a[0]*b[2]*c[3]-a[2]*b[3]*c[0]-a[3]*b[0]*c[2] \
	+a[0]*b[3]*c[2]+a[3]*b[2]*c[0]+a[2]*b[0]*c[3]
    y = -a[1]*b[0]*c[3]-a[0]*b[3]*c[1]-a[3]*b[1]*c[0] \
	+a[1]*b[3]*c[0]+a[3]*b[0]*c[1]+a[0]*b[1]*c[3]
    z = -a[1]*b[2]*c[0]-a[2]*b[0]*c[1]-a[0]*b[1]*c[2] \
	+a[1]*b[0]*c[2]+a[0]*b[2]*c[1]+a[2]*b[1]*c[0]
    return Vec4(t,-x,-y,-z)

def ComputePhi(pijt,pkt,pi):
    n_perp = pijt.Cross(pkt)
    if n_perp.P() == 0.0:
        n_perp = pijt.Cross(Vec4(0.,1.,0.,0.))
        if n_perp.P() == 0.0:
            n_perp = pijt.Cross(Vec4(0.,0.,1.,0.))
            if n_perp.P() == 0.0:
                n_perp = pijt.Cross(Vec4(0.,0.,0.,1.))
    n_perp /= n_perp.P()
    l_perp = LT(pijt,pkt,n_perp)
    l_perp /= m.sqrt(abs(l_perp.M2()))
    cp = -pi*n_perp
    sp = -pi*l_perp
    if sp==0 and cp==0: phi = 0.0
    else: phi = m.atan2(sp,cp)
    if cp<0.0:
        return phi+m.pi
    if sp>0.0:
        return phi
    return phi+2.0*m.pi

def ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk):
    pipj = pi*pj
    pipk = pi*pk
    pjpk = pj*pk
    yijk = pipj/(pipj+pipk+pjpk)
    zi = pipk/(pipk+pjpk)
    sij = (pi+pj).M2()
    Q2 = (pi+pj+pk).M2()
    po = pow(Q2-mij2-mk2,2)-4.0*mij2*mk2
    pn = pow(Q2-sij-mk2,2)-4.0*sij*mk2
    if pn<0.0 or po<0.0: return [0]
    Q = pi+pj+pk
    rpk = m.sqrt(po/pn)*(pk-(Q*pk)/Q2*Q)+(Q2+mk2-mij2)/(2.0*Q2)*Q
    rpi = Q-rpk
    res = [1,yijk,zi,ComputePhi(rpi,rpk,pi),rpi,rpk]
    return res

def ConstructFFDipole(mi2,mj2,mij2,mk2,pij,pk,ffp):
    n_perp = pij.Cross(pk)
    if n_perp.P() == 0.0:
        n_perp = pij.Cross(Vec4(0.,1.,0.,0.))
        if n_perp.P() == 0.0:
            n_perp = pij.Cross(Vec4(0.,0.,1.,0.))
            if n_perp.P() == 0.0:
                n_perp = pij.Cross(Vec4(0.,0.,0.,1.))
    n_perp /= n_perp.P();
    l_perp = LT(pij,pk,n_perp)
    l_perp /= m.sqrt(abs(l_perp.M2()))
    Q = pij+pk
    Q2 = Q.M2()
    sij = ffp[1]*(Q2-mk2)+(1.0-ffp[1])*(mi2+mj2)
    po = pow(Q2-mij2-mk2,2)-4.0*mij2*mk2
    pn = pow(Q2-sij-mk2,2)-4.0*sij*mk2
    if po<0.0 or pn<0.0: return [0]
    ecm = Q2-sij-mk2
    po = sign(ecm)*m.sqrt(po)
    pn = sign(ecm)*m.sqrt(pn)
    gam = 0.5*(ecm+pn)
    zt = ecm/pn*(ffp[2]-mk2/gam*(sij+mi2-mj2)/ecm)
    ktt = sij*zt*(1.0-zt)-(1.0-zt)*mi2-zt*mj2
    if ktt<0.0: return [0]
    ktt = m.sqrt(ktt)
    fpk = pn/po*(pk-(Q2-mij2+mk2)/(2.0*Q2)*Q)+(Q2-sij+mk2)/(2.0*Q2)*Q
    fpj = Q-fpk
    fpi = ktt*m.sin(ffp[3])*l_perp
    fpi += ktt*m.cos(ffp[3])*n_perp+zt/pn*(gam*fpj-sij*fpk)+ \
           (mi2+ktt*ktt)/zt/pn*(fpk-mk2/gam*fpj)
    fpj = Q-fpk-fpi
    if (fpi[0] > 0.) ^ (pij[0] > 0.): return [0]
    if (fpk[0] > 0.) ^ (pk[0] > 0.): return [0]
    if fpj[0] < 0: return [0]
    return [1,fpi,fpj,fpk]

def CorrectMomenta(p,mass):
    imax = 0
    emax = 0.0
    sum = -p[0]-p[1]
    for i in range(2,len(p)):
        sum+=p[i]
        if p[i][0]>emax:
            emax = p[i][0]
            imax = i
        p[i][0] = m.sqrt(p[i].P2()+mass[i]**2)
    p[imax] -= sum
    p[imax][0] = m.sqrt(p[imax].P2()+mass[imax]**2)
    e0tot = 0
    for i in range(len(p)):
        e0tot += -p[i][0] if i<2 else p[i][0]
    p2 = [ p[0].P2(), p[1].P2() ]
    E0 = [ -p[0][0], -p[1][0] ]
    E1 = [ p2[0]/E0[0], (p[0].Vec3()*p[1].Vec3())/E0[1] ]
    e1tot = E1[0]+E1[1]
    E2 = [ p2[0]*mass[0]**2/(2.*E0[0]**3), (p2[0]-E1[1]**2)/(2*E0[1]) ]
    e2tot = E2[0]+E2[1]
    eps1 = -e0tot/e1tot
    eps = eps1*(1-eps1*e2tot/e1tot)
    p[1] = p[1]-p[0]*eps
    p[0] = p[0]+p[0]*eps
    for i in range(2):
        if mass[i]==0.0:
            p[i][0]=p[i].P()
        else:
            p[i][0]=E0[i]+E1[i]*eps+E2[i]*eps**2

import sys, os, optparse, re
from xml.dom import minidom

parser = optparse.OptionParser()
parser.add_option("-f","--file",default="events.lhe",dest="infile")
parser.add_option("-o","--output",default="trajectories",dest="outpath")
parser.add_option("-i","--index-i",default=0,dest="idi")
parser.add_option("-j","--index-j",default=0,dest="idj")
parser.add_option("-k","--index-k",default=0,dest="idk")
parser.add_option("-t","--test",action="store_true",default=False,dest="check")
parser.add_option("-r","--range",default="[1e-3,1e-9,3]",dest="range")
parser.add_option("-p","--particles",default="[]",dest="pids")
parser.add_option("-c","--coll",action="store_true",default=False,dest="coll")
parser.add_option("-n","--numpoints",default=1000000,dest="npts")
parser.add_option("-C","--combine-is",action="store_true",default=False,dest="comb")
parser.add_option("-m","--masses",default="{}",dest="mass")
(opts,args) = parser.parse_args()

try:
    os.mkdir(opts.infile.split('.',1)[0])
except OSError:
    pass
idi = int(opts.idi)
idj = int(opts.idj)
idk = int(opts.idk)
if idj < 2: exit(1)
rang = eval(opts.range)
rang.append(1 if rang[2]<2 else rang[2]-1)
mass = eval(opts.mass)
pids = eval(opts.pids)
event = 0
events = minidom.parse(opts.infile)
points = file('{1}/{0}.dat'.format
              (opts.outpath,opts.infile.split('.',1)[0]),'w')
for moms in events.getElementsByTagName('event'):
    lines = moms.firstChild.nodeValue.split('\n')
    count = 0
    old = []
    ids = []
    for line in lines:
        line = re.sub(' +',' ',line).strip().split(' ')
        if line[0]=='' or len(line)<13: continue
        mom = Vec4(float(line[9]),float(line[6]),float(line[7]),float(line[8]))
        if count<2 and mom[0]>0.0: mom = -mom
        old.append(mom)
        ids.append(int(line[0]))
        count += 1
    for i, p in enumerate(old):
        if not i in mass: mass[i] = 0.
    event += 1
    if len(pids): ids = pids
    if opts.check:
        print 'point {0} {1}'.format(event,'{')
    Ecm2 = (old[0]+old[1]).M2()
    Q2 = (old[idi]+old[idj]+old[idk]).M2()
    for pt in range(rang[2]):
        lam = m.exp(m.log(rang[0])+pt/float(rang[3])*m.log(rang[1]/rang[0]))
        ffp = ClusterFFDipole(0.,0.,0.,0.,old[idi],old[idj],old[idk])
        if ffp[0]==0: continue
        pff = []
        for i in ffp: pff.append(i)
        if opts.coll:# coll
            pff[2] = pff[2]
            if idi > 1:
                if idk > 1: pff[1] = lam*Ecm2/Q2
                else: pff[1] = 1.0/(1.0+1.0/(lam*Ecm2/Q2/(1.0-pff[1])))
            else:
                if idk > 1: exit(1)
                else: pff[1] = -pff[2]*lam/(1.0-pff[2]*lam)
        else:# soft
            if (idi < 2) ^ (idk < 2): exit(1)
            C = pff[1]*pff[2]/(1.0-pff[2])
            if idi > 2:
                pff[2] = 1.0-2.0*lam/(1.0+C)
                pff[1] = C*(1.0-pff[2])/pff[2]
            else:
                pff[2] = 1.0/(1.0+2.0*(lam**2-lam*m.sqrt(1.0+lam**2))/(1.0+C))
                pff[1] = m.copysign(C*(1.0-pff[2])/pff[2],pff[1])
        cff = ConstructFFDipole(0.,0.,0.,0.,pff[4],pff[5],pff)
        if cff[0]==0: continue
        new = []
        for p in old: new.append(p)
        new[idi] = cff[1]
        new[idj] = cff[2]
        new[idk] = cff[3]
        newcms = -new[0]-new[1]
        zaxis = Rotation(newcms.Boost(new[0]),Vec4(1.,0.,0.,1.))
        for i, p in enumerate(new):
            new[i] = zaxis*newcms.Boost(p)
        new[0] = Vec4(-new[0][0],0.,0.,-new[0][3])
        new[1] = Vec4(-new[1][0],0.,0.,-new[1][3])
        CorrectMomenta(new,mass)
        if opts.check:
            new[0] = -new[0]
            new[1] = -new[1]
            cff = ClusterFFDipole(0.,0.,0.,0.,new[idi],new[idj],new[idk])
            new[0] = -new[0]
            new[1] = -new[1]
            sum = Vec4()
            for i, p in enumerate(new):
                sum += (p if i>(0 if opts.comb else 1) else -p)
                print "  id = {0}, p = {1}, p^2/p_0^2 = {2}".\
                    format(ids[i],p,p.M2()/p[0]**2)
            print "  y dev = {0}, z dev = {1}, mom dev = {2}".\
                format(cff[1]/pff[1]-1.,cff[2]/pff[2]-1.,sum)
        if opts.comb:
            new[1] = new[0]+new[1]
            new.pop(0)
        for i, p in enumerate(new):
            points.write("{0} {1:.16f} {2:.16f} {3:.16f} {4:.16f};\n".format
                         (ids[i],p[0],p[1],p[2],p[3]))
        points.write("\n")
    if opts.check==1:
        print '{0}'.format('}')
    if event == int(opts.npts): break
points.write("eof\n")
points.close()
