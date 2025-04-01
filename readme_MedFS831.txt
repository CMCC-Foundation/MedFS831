-- LOAD ENVIRONMENT --
module load anaconda/3-2022.10 impi-2021.6.0/2021.6.0 intel-2021.6.0/libszip/2.1.1-tvhyi oneapi-2022.1.0/2022.1.0 oneapi-2022.1.0/mkl/2022.1.0 intel-2021.6.0/impi-2021.6.0/hdf5-threadsafe/1.13.3-zbgha intel-2021.6.0/impi-2021.6.0/netcdf-c-threadsafe/4.9.0-wpe4t intel-2021.6.0/impi-2021.6.0/netcdf-fortran-threadsafe/4.6.0-75oow intel-2021.6.0/impi-2021.6.0/parallel-netcdf/1.12.3-eshb5 intel-2021.6.0/curl/7.85.0-djjip intel-2021.6.0/perl/5.36.0-jj4hw intel-2021.6.0/perl-uri/1.72-6at2i intel-2021.6.0/nco/5.0.6-jp6y4 intel-2021.6.0/cdo-threadsafe

install ttide_py from https://github.com/moflaher/ttide_py
install gsw-Python from https://teos-10.github.io/GSW-Python/
-- LOAD ENVIRONMENT --

01) git clone -b main git@github.com:CMCC-Foundation/MedFS831.git MedFS831 
02) cd MedFS831
03) git clone  --branch 4.2.0 https://forge.nemo-ocean.eu/nemo/nemo.git   medsea-nemo42
04) cp MyNEMO/MY_ARCH/arch-X64_JUNO.fcm medsea-nemo42/arch/.
05) -- EDIT: possibly modify arch-X64_JUNO.fcm according to computer architecture
06) cp MyNEMO/MY_INSTALL/run_makenemo.bsh medsea-nemo42/.
07) -- EDIT: XIOS path in medsea-nemo42/run_makenemo.bsh (vi medsea-nemo42/run_makenemo.bsh)
08) cd medsea-nemo42
09) sh run_makenemo.bsh
10) cp ../MyNEMO/MY_INSTALL/run_maketools.bsh tools/.
11) -- EDIT: XIOS path in tools/run_maketools.bsh (vi tools/run_maketools.bsh)
12) cd tools
13) sh run_maketools.bsh
14) cd DOMAINcfg
15) cp  ../../../MyNEMO/MY_NAMELIST/namelist_cfg_MedFS831_dom namelist_cfg
16) ln -fs ../../../DATA/STATIC/bathy_meter.nc .
17) ./make_domain_cfg.exe 
18) cd ../../../MyTOOLS/BDY_NEMO
19) python mk_bdycoords.py
20) cd ../../SCRIPT/
21) vi submit_MedFS831.bsh
21.1) date_end=20210101
21.2) OCEANVAR=false
21.3) ln_obsmisfits=false
21.4) bsub < submit_MedFS831.bsh
22) cd  ../MyTOOLS/COAST_DIST
23) ./compile.bsh
24) ./cstdst.x 
25) cd ../GRD_OCEANVAR
26) python mk_oceanvar_grd.py 
27) cd ../OBS_NEMO
28) bsub < run_prep.bsh
29) cd ../../
30) git clone -b main git@github.com:CMCC-Foundation/OceanVar2 oceanvar2 
31) cd  oceanvar2/src
32) ./link_in_work
33) cd ../work/
34) make all
35) cd ../../SCRIPT/
36) vi submit_MedFS831.bsh
36.1) date_end=20210131
36.2) OCEANVAR=true
36.3) ln_obsmisfits=true
36.4) bsub < submit_MedFS831.bsh 
37) vi submit_MedFS831.bsh
37.1) date_end=20210131
37.2) OCEANVAR=false
37.3) ln_obsmisfits=true
37.4) bsub < submit_MedFS831.bsh 
38) cd ../MyTOOLS/POST
39) vi submit_aggregate_NemoMisfits.bsh
39.1) analysis=false
39.2) bsub < submit_aggregate_NemoMisfits.bsh
40) vi submit_aggregate_NemoMisfits.bsh
40.1) analysis=true
40.2) bsub < submit_aggregate_NemoMisfits.bsh
41) python cmp_misfits_ins.py 
42) python cmp_misfits_sla.py 
43) python  mk_figure01.py
44) python  mk_figure04.py
45) python  mk_figure06.py

