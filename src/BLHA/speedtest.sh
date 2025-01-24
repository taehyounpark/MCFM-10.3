#!/bin/bash

nptsfac=1
while getopts n: flag
do
    case "${flag}" in
        n) nptsfac=${OPTARG};;
    esac
done

function do_test_ol() {
    olout=$(./oltest -e 10 -n $(echo $npts*$nptsfac | bc -l) $* | tail -n 15 | head -n 2 | sed -e "s/MCFM//" -e "s/OpenLoops//" -e "s/s//")
    rclout=0.0
    mlout=0.0
    echo $* $olout : $rclout : $mlout
}

function do_test_ol_rcl() {
    olout=$(./oltest -e 10 -n $(echo $npts*$nptsfac | bc -l) $* | tail -n 15 | head -n 2 | sed -e "s/MCFM//" -e "s/OpenLoops//" -e "s/s//")
    rclout=$(./rcltest -n $(echo $npts*$nptsfac | bc -l) $* | tail -n 14 | head -n 1 | sed -e "s/Recola//" -e "s/s//") #;s/(Ratio.*)//
    mlout=0.0
    echo $* $olout : $rclout : $mlout
}

function do_test_all() {
    olout=$(./oltest -e 10 -n $(echo $npts*$nptsfac | bc -l) $* | tail -n 15 | head -n 2 | sed -e "s/MCFM//" -e "s/OpenLoops//" -e "s/s//")
    rclout=$(./rcltest -n $(echo $npts*$nptsfac | bc -l) $* | tail -n 14 | head -n 1 | sed -e "s/Recola//" -e "s/s//") #;s/(Ratio.*)//
    mlout=$(./mltest -n $(echo $npts*$nptsfac | bc -l) $* | tail -n 14 | head -n 1 | sed -e "s/MadLoop//" -e "s/s//")
    echo $* $olout $rclout $mlout
}

function do_test_all_heft() {
    olout=$(./oltest -e 10 -n $(echo $npts*$nptsfac | bc -l) $* | tail -n 15 | head -n 2 | sed -e "s/MCFM//" -e "s/OpenLoops//" -e "s/s//")
    rclout=$(./rcltest_heft -n $(echo $npts*$nptsfac | bc -l) $* | tail -n 14 | head -n 1 | sed -e "s/Recola//" -e "s/s//") #;s/(Ratio.*)//
    mlout=$(./mltest -n $(echo $npts*$nptsfac | bc -l) $* | tail -n 14 | head -n 1 | sed -e "s/MadLoop//" -e "s/s//")
    echo $* $olout $rclout $mlout
}

function test_z() {
#echo -e "\033[1m===> Z processes <===\033[0m"
npts=10000
do_test_all d d~ e- e+
}
function test_zj() {
#echo -e "\033[1m===> Zj processes <===\033[0m"
npts=1000
do_test_all d d~ e- e+ g
}
function test_zjj() {
#echo -e "\033[1m===> Zjj processes <===\033[0m"
npts=100
do_test_all d d~ e- e+ g g
do_test_all d d~ e- e+ u u~
do_test_all d d~ e- e+ d d~
}
function test_wm() {
#echo -e "\033[1m===> W- processes <===\033[0m"
npts=10000
do_test_all u~ d e- ve~
}
function test_wmj() {
#echo -e "\033[1m===> W-j processes <===\033[0m"
npts=1000
do_test_all u~ d e- ve~ g
}
function test_wmjj() {
#echo -e "\033[1m===> W-jj processes <===\033[0m"
npts=100
do_test_all u~ d e- ve~ g g
do_test_all u~ d e- ve~ s s~
do_test_all u~ d e- ve~ u u~
do_test_all u~ d e- ve~ d d~
}
function test_tt() {
#echo -e "\033[1m===> tt processes <===\033[0m"
npts=10000
do_test_all d d~ t t~
do_test_all g g t t~
}
function test_heft() {
#echo -e "\033[1m===> h (EFT) processes <===\033[0m"
npts=100000
do_test_all_heft g g h -Pmodel=heft
}
function test_hjeft() {
#echo -e "\033[1m===> hj (EFT) processes <===\033[0m"
npts=10000
do_test_all_heft g g h g -Pmodel=heft
do_test_all_heft d d~ h g -Pmodel=heft
}
function test_hjjeft() {
#echo -e "\033[1m===> hjj (EFT) processes <===\033[0m"
npts=1000
do_test_all_heft g g h g g -Pmodel=heft
do_test_all_heft d d~ h g g -Pmodel=heft
do_test_all_heft d d~ h d d~ -Pmodel=heft
do_test_all_heft d d~ h u u~ -Pmodel=heft
}
function test_hsm() {
#echo -e "\033[1m===> h (SM) processes <===\033[0m"
npts=100000
do_test_all g g h
}
function test_hjsm() {
#echo -e "\033[1m===> hj (SM) processes <===\033[0m"
npts=10000
do_test_all g g h g
do_test_all d d~ h g
}
function test_hjjsm() {
#echo -e "\033[1m===> hjj (SM) processes <===\033[0m"
npts=1000
do_test_all g g h g g
do_test_all d d~ h g g
do_test_all d d~ h d d~
do_test_all d d~ h u u~
}
function test_hhsm() {
#echo -e "\033[1m===> hh (SM) processes <===\033[0m"
npts=10000
do_test_all g g h h -S2 -W2 -A12
}
function test_yj() {
#echo -e "\033[1m===> yj processes <===\033[0m"
npts=10000
do_test_all d d~ a g
}
function test_yjj() {
#echo -e "\033[1m===> yjj processes <===\033[0m"
npts=1000
do_test_all d~ d a g g -Ptop_mass=100000 -Ptop_yukawa=100000
do_test_all d d~ a d d~ -Ptop_mass=100000 -Ptop_yukawa=100000
do_test_all d d~ a u u~ -Ptop_mass=100000 -Ptop_yukawa=100000
}
function test_yy() {
#echo -e "\033[1m===> yy processes <===\033[0m"
npts=10000
do_test_all d d~ a a
do_test_all g g a a -S2 -W2 -A12
}
function test_yyj() {
#echo -e "\033[1m===> yyj processes <===\033[0m"
npts=1000
do_test_all d d~ a a g -Ptop_mass=100000 -Ptop_yukawa=100000
}
function test_yyy() {
#echo -e "\033[1m===> yyy processes <===\033[0m"
npts=100
do_test_all d d~ a a a
}
function test_yyyy() {
#echo -e "\033[1m===> yyy processes <===\033[0m"
npts=100
do_test_all d d~ a a a a
}
function test_ww() {
#echo -e "\033[1m===> WW processes <===\033[0m"
npts=1000
do_test_all d d~ e- ve~ mu+ vmu
}
function test_wz() {
#echo -e "\033[1m===> WZ processes <===\033[0m"
npts=1000
do_test_all d~ u e- e+ mu+ vmu
do_test_all d u~ e- ve~ vmu vmu~
}
function test_zz() {
#echo -e "\033[1m===> ZZ processes <===\033[0m"
npts=1000
do_test_all d d~ e- e+ mu- mu+
do_test_all d d~ e- e+ vmu vmu~
do_test_all d d~ ve ve~ mu+ mu-
}
function test_wy() {
#echo -e "\033[1m===> Wy processes <===\033[0m"
npts=1000
do_test_all d u~ e- ve~ a
}
function test_zy() {
#echo -e "\033[1m===> Zy processes <===\033[0m"
npts=1000
do_test_all d d~ e- e+ a
do_test_all d d~ ve ve~ a
}
function test_wh() {
#echo -e "\033[1m===> Wh processes <===\033[0m"
npts=10000
do_test_all d u~ e- ve~ h
}
function test_whj() {
#echo -e "\033[1m===> Whj processes <===\033[0m"
npts=1000
do_test_all d u~ e- ve~ h g
}
function test_zh() {
#echo -e "\033[1m===> Zh processes <===\033[0m"
npts=10000
do_test_all d d~ e- e+ h
}
function test_zhj() {
#echo -e "\033[1m===> Zhj processes <===\033[0m"
npts=1000
do_test_all d d~ e- e+ h g
}
function test_jj() {
#echo -e "\033[1m===> jj processes <===\033[0m"
npts=10000
do_test_all g g g g -Ptop_mass=100000 -Ptop_yukawa=100000
do_test_all d d~ d d~
do_test_all d d~ u u~
do_test_all d d~ g g 
}


grep -m 1 'model name' /proc/cpuinfo | sed -e "s/CPU//" -e "s/model name//" -e "s/://" -e "s/(R)//g"
test_z
test_zj
test_zjj
test_wm
test_wmj
test_wmjj
test_wy
test_zy
test_ww
test_wz
test_zz
test_wh
test_whj
test_zh
test_zhj
test_yj
test_yjj
test_yy
test_yyj
test_yyy
test_yyyy
test_jj
test_tt
test_hsm
test_hjsm
test_hjjsm
test_hhsm
test_heft
test_hjeft
test_hjjeft


