

c----These are reduced coefficients
      double complex coeff(0:4,60),scints(4,60,-2:0)


      integer,parameter::d234x56=1,d27x56x34=2,d34x156=3,d34x56x27=4,
     & d127=5,
     & d34x56x27sl=6,d3456sla=7,d3456slb=8,d56x34x17sl=9,
     & d56x234sla=10,d56x234slb=11,d34x156sla=12,d34x156slb=13,
     & d34x567sl=14,d34x12x56sla=15,d34x12x56slb=16,
     & d56x347sl=17,d12x56x34sl=18,
     & d7x1x56=1,d1x27x34=2,d7x2x34=3,d2x17x56=4,d1x7x2=5,
     & d23x1x45q=6,d23x1x67q=7,d45x1x67q=8,d17x2x56=19,d27x1x34=20,
     & d7x1x34=21,d7x2x56=22,d12x7x34=23,d1x2x56=24,d1x34x7=25,
     & d1x7x34=26,d2x1x34=27,d2x1x7=28,d2x56x7=29,d2x7x56=30,
     & d1x7x34sl=31,d12x7x34sl=32,d17x2x56sl=33,d1x2x56sl=34,
     & d1x34x7sl=35,d27x1x34sl=36,d2x1x34sl=37,d2x56x7sl=38,d2x7x56sl=39


      integer,parameter::c56x17=1,c34x27=2,c56x34=3,c12x34sl=4,
     & c12x56sl=5,c56x17sl=6,c34x27sl=7,c56x34sl=8
      integer,parameter::
     &  c12x34=9,c1x27=10,c1x56=11,c1x7=12,c2x17=13,c2x34=14,c2x7=15,
     &  c1x27sl=17,
     &  c1x2sl=18,
     &  c1x567sl=19,
     &  c1x56sl=20,
     &  c1x7sl=21,
     &  c2x17sl=22,
     &  c2x347sl=23,
     &  c2x34sl=24,
     &  c7x12sl=27,
     &  c7x234sl=28,
     &  c7x34sl=29,
     &  c7x56sl=30,

     & c1x23q=31,
     & c1x45q=32,
     & c1x67q=33,
     & c23x45q=34,
     & c23x67q=35,
     & c45x67q=36,
     & c17x34=37,c1x34=38,c27x56=39,c2x56=40,c7x134=41,c17x34sl=42,
     & c27x56sl=43,
     & c17x56sl=6,
     &  c1x347sl=45,
     &  c1x34sl=46,
     &  c2x134sl=47,
     &  c2x56sl=48,
     &  c7x156=49,
     &  c34x56sl=8,
     &  c7x134sl=50

      integer,parameter::
     & b127=1,b17=2,b27=3,b156=4,b234=5,b56=6,b34=7,b256=8,b134=9,
     & b17sl=10,b12sl=11,b127sl=12,b567sl=13,b347sl=14,
     & b56sl=15,b34sl=16,b27sl=17,
     & b156sl=18,b256sl=19,b134sl=20,b234sl=21,
     & b123q=22,b145q=23,b167q=24,b23q=25,b45q=26,b67q=27

      integer irat,iratsl
      parameter (irat=1,iratsl=2)

      integer,parameter::
c d1x27x34sl.f:!     d34x56x27sl=6, formerly box3sl
     & d1x27x34sl=6,
c d2x1x7sl.f:!     d3456sla=7, formerly box8sl
     & d2x1x7sl=7,
c d1x2x7sl.f:!     d3456slb=8, formerly box7sl
     & d1x2x7sl=8,
c d56x34x17sl.f:!      d56x34x17sl=9, formerly box6sl
     & d17x2x34sl=9,
c d1x56x7sl.f:!     d56x234sla=10,formerly box1sl
     & d1x56x7sl=10,
c d1x7x56sl.f:!     d56x234slb=11,formerly box2sl
     & d1x7x56sl=11,
c d2x34x7sl.f:!     d34x156sla=12, formerly box4sl
     & d2x34x7sl=12,
c d2x7x34sl.f:!     d34x156slb=13, i.e. formerly box5sl
     & d2x7x34sl=13,
c d1x2x34sl.f:!      d34x567sl=14, formerly box9sl
     & d1x2x34sl=14,
c d7x34x12sl.f:!     d34x12x56sla=15, formerly box10sl
     & d7x34x12sl=15,
c d2x1x56sl.f:!     This is d56x347sl=17 formerly box11sl
     & d2x1x56sl=17,
c d7x12x34sl.f:!     This is known as d34x12x56slb=16 formerly box12sl
     & d7x12x34sl=16,
c d7x12x56sl.f:!     d12x56x34sl=18, formerly box13sl
     & d7x12x56sl=18

