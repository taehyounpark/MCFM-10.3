target_sources(objlib PRIVATE
AP0.f
AP1.f
APP0xP0.f
APtI1xP0.f
calI1.f
calI2.f
CMakeLists.txt
Gammafill.f
hard_DY.f
hard_gamgam.f
hard_H.f
hard_Vgam.f
I2slash.f
qtassemble.f
qtbeam1.f
qtbeam2.f
qtint.f
qtlumxmsq.f
qtsoft.f
qt_xmsq.f
tildeI1.f
tildeI2.f
tildep0.f
tildeSb1.f
tildeSb2.f
)

if(with_vvamp)
    target_sources(objlib PRIVATE
        hard_VV.f
    )
endif(with_vvamp)
