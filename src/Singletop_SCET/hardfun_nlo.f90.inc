    w(1)=Hc1(0)
    w(2)=Hc1(1)
    w(3)=Hc2(1,0)
    w(4)=1/(1 - omx)
    w(5)=1/(mt - mt*omx)
    w(6)=2*w(2)
    w(7)= - w(6) + w(4) - 3 - 2*w(1)
    w(7)=w(7)*w(6)
    w(8)=2*LogMtMu
    w(9)= - w(8) + 5 + 4*w(2)
    w(8)=w(9)*w(8)
    w(9)=Pi**2
    w(7)=w(8) + w(7) + 4*w(3) - 12 - 1.D0/6.D0*w(9)
    w(7)=tree*w(7)
    w(6)=w(6)*w(5)
    w(8)= - struc*w(6)
    w(7)=w(8) + w(7)
    w(7)=treeCC*w(7)
    w(6)= - tree*strucCC*w(6)
    w(6)=w(6) + w(7)
    w(6)=CF*w(6)
    w(7)=treeCC*tree


    hardfun(0) = w(7)

    hardfun(1) = w(6)

    hardfun(2) =  0
