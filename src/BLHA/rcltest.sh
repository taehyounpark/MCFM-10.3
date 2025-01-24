#!/bin/bash

target=$1
test -z "$target" && target=6
echo "Target accuracy: "$target" digits"

function do_test(){
  if (( $target == 0 )); then
    res=$(./rcltest -pn1 $* | grep MCFM\ = | sed 's/MCFM.*Ratio =//1;s/Finite:/F =/1;s/IR:/IR =/;s/IR2:/IR2 =/;s/Born:/B =/')
  else
    res=$(./rcltest -pn1 $* | grep MCFM\ = | sed 's/MCFM.*Ratio =//1;s/Finite:/F =/1;s/IR:/IR =/;s/IR2:/IR2 =/;s/Born:/B =/' |\
    awk '{for(i=1;i<=NF;i+=1){a[2]=a[1];a[1]=a[0];a[0]=$i;if(a[1]=="="){\
          dev=a[0]-1;if(dev>-1e-'$target'&&dev<1e-'$target')print a[2]" \033[32mok\033[0m";else print "\033[31mfailed\033[0m";}}}')
  fi
  echo $*": "$res
}
function do_test_heft(){
  if (( $target == 0 )); then
    res=$(./rcltest_heft -pn1 $* | grep MCFM\ = | sed 's/MCFM.*Ratio =//1;s/Finite:/F =/1;s/IR:/IR =/;s/IR2:/IR2 =/;s/Born:/B =/')
  else
    res=$(./rcltest_heft -pn1 $* | grep MCFM\ = | sed 's/MCFM.*Ratio =//1;s/Finite:/F =/1;s/IR:/IR =/;s/IR2:/IR2 =/;s/Born:/B =/' |\
    awk '{for(i=1;i<=NF;i+=1){a[2]=a[1];a[1]=a[0];a[0]=$i;if(a[1]=="="){\
          dev=a[0]-1;if(dev>-1e-'$target'&&dev<1e-'$target')print a[2]" \033[32mok\033[0m";else print "\033[31mfailed\033[0m";}}}')
  fi
  echo $*": "$res
}

function test_z() {
echo -e "\033[1m===> Z processes <===\033[0m"
do_test d d~ e- e+
do_test d~ d e- e+
do_test u u~ e- e+
do_test u~ u e- e+
}
function test_zj() {
echo -e "\033[1m===> Zj processes <===\033[0m"
do_test d d~ e- e+ g
do_test d g e- e+ d
do_test d~ d e- e+ g
do_test d~ g e- e+ d~
do_test u u~ e- e+ g
do_test u g e- e+ u
do_test u~ u e- e+ g
do_test u~ g e- e+ u~
do_test g d e- e+ d
do_test g d~ e- e+ d~
do_test g u e- e+ u
do_test g u~ e- e+ u~
}
function test_zjj() {
echo -e "\033[1m===> Zjj processes <===\033[0m"
do_test d d e- e+ d d
do_test d d~ e- e+ d d~
do_test d d~ e- e+ u u~
do_test d d~ e- e+ s s~
do_test d d~ e- e+ g g
do_test d u e- e+ d u
do_test d u~ e- e+ d u~
do_test d s e- e+ d s
do_test d s~ e- e+ d s~
do_test d g e- e+ g d
do_test d~ d e- e+ d d~
do_test d~ d e- e+ u u~
do_test d~ d e- e+ s s~
do_test d~ d e- e+ g g
do_test d~ d~ e- e+ d~ d~
do_test d~ u e- e+ u d~
do_test d~ u~ e- e+ d~ u~
do_test d~ s e- e+ s d~
do_test d~ s~ e- e+ d~ s~
do_test d~ g e- e+ g d~
do_test u d e- e+ d u
do_test u d~ e- e+ u d~
do_test u u e- e+ u u
do_test u u~ e- e+ d d~
do_test u u~ e- e+ u u~
do_test u u~ e- e+ c c~
do_test u u~ e- e+ g g
do_test u s e- e+ u s
do_test u c e- e+ u c
do_test u c~ e- e+ u c~
do_test u g e- e+ g u
do_test u~ d e- e+ d u~
do_test u~ d~ e- e+ d~ u~
do_test u~ u e- e+ d d~
do_test u~ u e- e+ u u~
do_test u~ u e- e+ c c~
do_test u~ u e- e+ g g
do_test u~ u~ e- e+ u~ u~
do_test u~ s~ e- e+ u~ s~
do_test u~ c e- e+ c u~
do_test u~ c~ e- e+ u~ c~
do_test u~ g e- e+ g u~
do_test s d e- e+ d s
do_test s u e- e+ u s
do_test s~ d~ e- e+ d~ s~
do_test s~ u~ e- e+ u~ s~
do_test c u e- e+ u c
do_test c~ u~ e- e+ u~ c~
do_test g d e- e+ g d
do_test g d~ e- e+ g d~
do_test g u e- e+ g u
do_test g u~ e- e+ g u~
do_test g g e- e+ d d~
do_test g g e- e+ u u~
}
function test_wm() {
echo -e "\033[1m===> W- processes <===\033[0m"
do_test d u~ e- ve~
do_test u~ d e- ve~
}
function test_wmj() {
echo -e "\033[1m===> W-j processes <===\033[0m"
do_test d u~ e- ve~ g
do_test d g e- ve~ u
do_test u~ d e- ve~ g
do_test u~ g e- ve~ d~
do_test g d e- ve~ u
do_test g u~ e- ve~ d~
}
function test_wmjj() {
echo -e "\033[1m===> W-jj processes <===\033[0m"
do_test d d e- ve~ d u
do_test d d~ e- ve~ u d~
do_test d d~ e- ve~ c s~
do_test d u e- ve~ u u
do_test d u~ e- ve~ d d~
do_test d u~ e- ve~ u u~
do_test d u~ e- ve~ s s~
do_test d u~ e- ve~ g g
do_test d s e- ve~ d c
do_test d s e- ve~ u s
do_test d s~ e- ve~ u s~
do_test d c~ e- ve~ d s~
do_test d g e- ve~ g u
do_test d~ d e- ve~ u d~
do_test d~ d e- ve~ c s~
do_test d~ u~ e- ve~ d~ d~
do_test d~ s e- ve~ c d~
do_test d~ c~ e- ve~ d~ s~
do_test u d e- ve~ u u
do_test u u~ e- ve~ u d~
do_test u~ d e- ve~ d d~
do_test u~ d e- ve~ u u~
do_test u~ d e- ve~ s s~
do_test u~ d e- ve~ g g
do_test u~ d~ e- ve~ d~ d~
do_test u~ u e- ve~ u d~
do_test u~ u~ e- ve~ d~ u~
do_test u~ s e- ve~ s d~
do_test u~ s~ e- ve~ d~ s~
do_test u~ g e- ve~ g d~
do_test s d e- ve~ d c
do_test s d e- ve~ u s
do_test s~ u~ e- ve~ d~ s~
do_test c~ d~ e- ve~ d~ s~
do_test g d e- ve~ g u
do_test g u~ e- ve~ g d~
do_test g g e- ve~ u d~
}
function test_wp() {
echo -e "\033[1m===> W+ processes <===\033[0m"
do_test d~ u e+ ve
do_test u d~ e+ ve
}
function test_wpj() {
echo -e "\033[1m===> W+j processes <===\033[0m"
do_test d~ u e+ ve g
do_test d~ g e+ ve u~
do_test u d~ e+ ve g
do_test u g e+ ve d
do_test g d~ e+ ve u~
do_test g u e+ ve d
}
function test_wpjj() {
echo -e "\033[1m===> W+jj processes <===\033[0m"
do_test d d~ e+ ve d u~
do_test d d~ e+ ve s c~
do_test d u e+ ve d d
do_test d s~ e+ ve d c~
do_test d c e+ ve d s
do_test d~ d e+ ve d u~
do_test d~ d e+ ve s c~
do_test d~ d~ e+ ve d~ u~
do_test d~ u e+ ve d d~
do_test d~ u e+ ve u u~
do_test d~ u e+ ve s s~
do_test d~ u e+ ve g g
do_test d~ u~ e+ ve u~ u~
do_test d~ s e+ ve s u~
do_test d~ s~ e+ ve d~ c~
do_test d~ s~ e+ ve u~ s~
do_test d~ c e+ ve s d~
do_test d~ g e+ ve g u~
do_test u d e+ ve d d
do_test u d~ e+ ve d d~
do_test u d~ e+ ve u u~
do_test u d~ e+ ve s s~
do_test u d~ e+ ve g g
do_test u u e+ ve d u
do_test u u~ e+ ve d u~
do_test u s e+ ve d s
do_test u s~ e+ ve d s~
do_test u g e+ ve g d
do_test u~ d~ e+ ve u~ u~
do_test u~ u e+ ve d u~
do_test s u e+ ve d s
do_test s~ d~ e+ ve d~ c~
do_test s~ d~ e+ ve u~ s~
do_test c d e+ ve d s
do_test g d~ e+ ve g u~
do_test g u e+ ve g d
do_test g g e+ ve d u~
}
function test_tt() {
echo -e "\033[1m===> tt processes <===\033[0m"
do_test d d~ t t~
do_test d~ d t t~
do_test g g t t~
}
function test_heft() {
echo -e "\033[1m===> h (EFT) processes <===\033[0m"
do_test_heft g g h -Pmodel=heft
}
function test_hjeft() {
echo -e "\033[1m===> hj (EFT) processes <===\033[0m"
do_test_heft d d~ h g -Pmodel=heft
do_test_heft d g h d -Pmodel=heft
do_test_heft d~ d h g -Pmodel=heft
do_test_heft d~ g h d~ -Pmodel=heft
do_test_heft g d h d -Pmodel=heft
do_test_heft g d~ h d~ -Pmodel=heft
do_test_heft g g h g -Pmodel=heft
}
function test_hjjeft() {
echo -e "\033[1m===> hjj (EFT) processes <===\033[0m"
do_test_heft d d h d d -Pmodel=heft
do_test_heft d d~ h d d~ -Pmodel=heft
do_test_heft d d~ h u u~ -Pmodel=heft
do_test_heft d d~ h g g -Pmodel=heft
do_test_heft d u h d u -Pmodel=heft
do_test_heft d u~ h d u~ -Pmodel=heft
do_test_heft d g h g d -Pmodel=heft
do_test_heft d~ d h d d~ -Pmodel=heft
do_test_heft d~ d h u u~ -Pmodel=heft
do_test_heft d~ d h g g -Pmodel=heft
do_test_heft d~ d~ h d~ d~ -Pmodel=heft
do_test_heft d~ u h u d~ -Pmodel=heft
do_test_heft d~ u~ h d~ u~ -Pmodel=heft
do_test_heft d~ g h g d~ -Pmodel=heft
do_test_heft u d h d u -Pmodel=heft
do_test_heft u~ d~ h d~ u~ -Pmodel=heft
do_test_heft g d h g d -Pmodel=heft
do_test_heft g d~ h g d~ -Pmodel=heft
do_test_heft g g h d d~ -Pmodel=heft
do_test_heft g g h g g -Pmodel=heft
}
function test_hsm() {
echo -e "\033[1m===> h (SM) processes <===\033[0m"
do_test g g h
}
function test_hjsm() {
echo -e "\033[1m===> hj (SM) processes <===\033[0m"
do_test d d~ h g
do_test d g h d
do_test d~ d h g
do_test d~ g h d~
do_test g d h d
do_test g d~ h d~
do_test g g h g
}
function test_hjjsm() {
echo -e "\033[1m===> hjj (SM) processes <===\033[0m"
do_test d d h d d
do_test d d~ h d d~
do_test d d~ h u u~
do_test d d~ h g g
do_test d u h d u
do_test d u~ h d u~
do_test d g h g d
do_test d~ d h d d~
do_test d~ d h u u~
do_test d~ d h g g
do_test d~ d~ h d~ d~
do_test d~ u h u d~
do_test d~ u~ h d~ u~
do_test d~ g h g d~
do_test u d h d u
do_test u~ d~ h d~ u~
do_test g d h g d
do_test g d~ h g d~
do_test g g h d d~
do_test g g h g g
}
function test_yj() {
echo -e "\033[1m===> yj processes <===\033[0m"
do_test d d~ a g
do_test d g a d
do_test d~ d a g
do_test d~ g a d~
do_test u u~ a g
do_test u g a u
do_test u~ u a g
do_test u~ g a u~
do_test g d a d
do_test g d~ a d~
do_test g u a u
do_test g u~ a u~
}
function test_yjj() {
echo -e "\033[1m===> yjj processes <===\033[0m"
do_test d d a d d -Ptop_mass=100000
do_test d d~ a d d~ -Ptop_mass=100000
do_test d d~ a u u~ -Ptop_mass=100000
do_test d d~ a s s~ -Ptop_mass=100000
do_test d d~ a g g -Ptop_mass=100000
do_test d u a d u -Ptop_mass=100000
do_test d u~ a d u~ -Ptop_mass=100000
do_test d s a d s -Ptop_mass=100000
do_test d s~ a d s~ -Ptop_mass=100000
do_test d g a g d -Ptop_mass=100000
do_test d~ d a d d~ -Ptop_mass=100000
do_test d~ d a u u~ -Ptop_mass=100000
do_test d~ d a s s~ -Ptop_mass=100000
do_test d~ d a g g -Ptop_mass=100000
do_test d~ d~ a d~ d~ -Ptop_mass=100000
do_test d~ u a u d~ -Ptop_mass=100000
do_test d~ u~ a d~ u~ -Ptop_mass=100000
do_test d~ s a s d~ -Ptop_mass=100000
do_test d~ s~ a d~ s~ -Ptop_mass=100000
do_test d~ g a g d~ -Ptop_mass=100000
do_test u d a d u -Ptop_mass=100000
do_test u d~ a u d~ -Ptop_mass=100000
do_test u u a u u -Ptop_mass=100000
do_test u u~ a d d~ -Ptop_mass=100000
do_test u u~ a u u~ -Ptop_mass=100000
do_test u u~ a c c~ -Ptop_mass=100000
do_test u u~ a g g -Ptop_mass=100000
do_test u s a u s -Ptop_mass=100000
do_test u c a u c -Ptop_mass=100000
do_test u c~ a u c~ -Ptop_mass=100000
do_test u g a g u -Ptop_mass=100000
do_test u~ d a d u~ -Ptop_mass=100000
do_test u~ d~ a d~ u~ -Ptop_mass=100000
do_test u~ u a d d~ -Ptop_mass=100000
do_test u~ u a u u~ -Ptop_mass=100000
do_test u~ u a c c~ -Ptop_mass=100000
do_test u~ u a g g -Ptop_mass=100000
do_test u~ u~ a u~ u~ -Ptop_mass=100000
do_test u~ s~ a u~ s~ -Ptop_mass=100000
do_test u~ c a c u~ -Ptop_mass=100000
do_test u~ c~ a u~ c~ -Ptop_mass=100000
do_test u~ g a g u~ -Ptop_mass=100000
do_test s d a d s -Ptop_mass=100000
do_test s u a u s -Ptop_mass=100000
do_test s~ d~ a d~ s~ -Ptop_mass=100000
do_test s~ u~ a u~ s~ -Ptop_mass=100000
do_test c u a u c -Ptop_mass=100000
do_test c~ u~ a u~ c~ -Ptop_mass=100000
do_test g d a g d -Ptop_mass=100000
do_test g d~ a g d~ -Ptop_mass=100000
do_test g u a g u -Ptop_mass=100000
do_test g u~ a g u~ -Ptop_mass=100000
do_test g g a d d~ -Ptop_mass=100000
do_test g g a u u~ -Ptop_mass=100000
}
function test_yy() {
echo -e "\033[1m===> yy processes <===\033[0m"
do_test d d~ a a
do_test d~ d a a
do_test u u~ a a
do_test u~ u a a
do_test g g a a -S2 -W2 -A12
}
function test_yyj() {
echo -e "\033[1m===> yyj processes <===\033[0m"
do_test d d~ a a g -Ptop_mass=100000
do_test d g a a d -Ptop_mass=100000
do_test d~ d a a g -Ptop_mass=100000
do_test d~ g a a d~ -Ptop_mass=100000
do_test u u~ a a g -Ptop_mass=100000
do_test u g a a u -Ptop_mass=100000
do_test u~ u a a g -Ptop_mass=100000
do_test u~ g a a u~ -Ptop_mass=100000
do_test g d a a d -Ptop_mass=100000
do_test g d~ a a d~ -Ptop_mass=100000
do_test g u a a u -Ptop_mass=100000
do_test g u~ a a u~ -Ptop_mass=100000
}
function test_yyy() {
echo -e "\033[1m===> yyy processes <===\033[0m"
do_test d d~ a a a
do_test d~ d a a a
do_test u u~ a a a
do_test u~ u a a a
}
function test_yyyy() {
echo -e "\033[1m===> yyyy processes <===\033[0m"
do_test d d~ a a a a
do_test d~ d a a a a
do_test u u~ a a a a
do_test u~ u a a a a
}
function test_ww() {
echo -e "\033[1m===> WW processes <===\033[0m"
do_test d d~ e- ve~ mu+ vmu
do_test u u~ e- ve~ mu+ vmu
do_test d d~ e+ ve mu- vmu~
do_test u u~ e+ ve mu- vmu~
}
function test_wz() {
echo -e "\033[1m===> WZ processes <===\033[0m"
do_test d~ u e- e+ mu+ vmu
do_test u d~ e- e+ mu+ vmu
do_test d u~ e- mu- e+ vmu~
do_test u~ d e- mu- e+ vmu~
do_test d u~ e- vmu ve~ vmu~
do_test u~ d e- vmu ve~ vmu~
do_test d~ u e+ ve vmu vmu~
do_test u d~ e+ ve vmu vmu~
}
function test_zz() {
echo -e "\033[1m===> ZZ processes <===\033[0m"
do_test d d~ e- mu- e+ mu+
do_test d d~ e- e- e+ e+
do_test d d~ e- e+ vmu vmu~
do_test d d~ ve ve~ mu+ mu-
}
function test_wy() {
echo -e "\033[1m===> Wy processes <===\033[0m"
do_test d u~ e- ve~ a
do_test u~ d e- ve~ a
do_test d~ u e+ ve a
do_test u d~ e+ ve a
}
function test_zy() {
echo -e "\033[1m===> Zy processes <===\033[0m"
do_test d d~ e- e+ a
do_test d~ d e- e+ a
do_test u u~ e- e+ a
do_test u~ u e- e+ a
do_test d d~ ve ve~ a
do_test d~ d ve ve~ a
do_test u u~ ve ve~ a
do_test u~ u ve ve~ a
}
function test_wh() {
echo -e "\033[1m===> Wh processes <===\033[0m"
do_test d u~ e- ve~ h
do_test u~ d e- ve~ h
do_test d~ u e+ ve h
do_test u d~ e+ ve h
}
function test_whj() {
echo -e "\033[1m===> Whj processes <===\033[0m"
do_test d u~ e- ve~ h g
do_test d g e- ve~ h u
do_test u~ d e- ve~ h g
do_test u~ g e- ve~ h d~
do_test g d e- ve~ h u
do_test g u~ e- ve~ h d~
do_test d~ u e+ ve h g
do_test d~ g e+ ve h u~
do_test u d~ e+ ve h g
do_test u g e+ ve h d
do_test g d~ e+ ve h u~
do_test g u e+ ve h d
}
function test_zh() {
echo -e "\033[1m===> Zh processes <===\033[0m"
do_test d d~ e- e+ h
do_test d~ d e- e+ h
do_test u u~ e- e+ h
do_test u~ u e- e+ h
}
function test_zhj() {
echo -e "\033[1m===> Zhj processes <===\033[0m"
do_test d d~ e- e+ h g -Ptop_mass=100000
do_test d g e- e+ h d -Ptop_mass=100000
do_test d~ d e- e+ h g -Ptop_mass=100000
do_test d~ g e- e+ h d~ -Ptop_mass=100000
do_test u u~ e- e+ h g -Ptop_mass=100000
do_test u g e- e+ h u -Ptop_mass=100000
do_test u~ u e- e+ h g -Ptop_mass=100000
do_test u~ g e- e+ h u~ -Ptop_mass=100000
do_test g d e- e+ h d -Ptop_mass=100000
do_test g d~ e- e+ h d~ -Ptop_mass=100000
do_test g u e- e+ h u -Ptop_mass=100000
do_test g u~ e- e+ h u~ -Ptop_mass=100000
}
function test_jj() {
echo -e "\033[1m===> jj processes <===\033[0m"
do_test d d d d
do_test d d~ d d~
do_test d d~ u u~
do_test d d~ g g
do_test d u d u
do_test d u~ d u~
do_test d g g d
do_test d~ d d d~
do_test d~ d u u~
do_test d~ d g g
do_test d~ d~ d~ d~
do_test d~ u u d~
do_test d~ u~ d~ u~
do_test d~ g g d~
do_test u d d u
do_test u~ d~ d~ u~
do_test g d g d
do_test g d~ g d~
do_test g g d d~
do_test g g g g -Ptop_mass=100000
}

test_z
test_zj
test_zjj
test_wm
test_wmj
test_wmjj
test_wp
test_wpj
test_wpjj
test_tt
test_heft
test_hjeft
test_hjjeft
test_hsm
test_hjsm
test_hjjsm
test_yj
test_yjj
test_yy
test_yyj
test_yyy
test_yyyy
test_wy
test_zy
test_ww
test_wz
test_zz
test_wh
test_whj
test_zh
test_zhj
test_jj
