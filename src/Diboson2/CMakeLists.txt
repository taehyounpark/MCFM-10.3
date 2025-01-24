target_sources(objlib PRIVATE
fillcoupfac.f
Li22.f
MXXX.f
Qformfill.f
qqbAjF.f
qqbVVMLLL.f
qqbVVMRLL.f
qqbVVMXXXtreeS.f
qt_amp.f
VVampfill.f
vv_msq.f
)

if(with_vvamp)
    target_sources(objlib PRIVATE
        vv_amps.f
        lumxmsq_vv.f
    )
endif(with_vvamp)
