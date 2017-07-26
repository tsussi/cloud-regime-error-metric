region=['2020','exts','extn','icey']
refrgn=['2020','extra','extra','icey']
nregimes=[7,7,7,6]
; model data info
model_names=['MOHC','CCCma','CNRM-CERFACS','MRI','MIROC']
models=['HadGEM2-A','CanAM4','CNRM-CM5','MRI-CGCM3','MIROC5']
versions=['v20111105','v2','v20120203','v1','v20120201']
period1=['19780901','19790101','19790101','19790101','19790101']
period2=['20081230','19791231','20081231','20100228','20081231']

nmodels=n_elements(models)
print, 'nmodels',nmodels

nregions=n_elements(region)
print, 'nregions',nregions

; model data directory and file
; HadGEM2 does not provide isicemask data
;for m=0,nmodels-1 do begin
for m=1,1 do begin
tgsic=ncassoc('/data/local2/hadyt/cmip5/CanAM4/sic_day_CanAM4_amip_r1i1p1_19790101-19791231.nc',/ignore)
tgsnc=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/snc_day_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
tgalbis=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/albisccp_cfDay_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
tgpctis=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/pctisccp_cfDay_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
tgcltis=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/cltisccp_cfDay_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
tgrsdt=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/rsdt_cfDay_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
tgnrsdt=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/rsdt_cfDay_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
tgtsw=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/rsut_cfDay_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
tgtlw=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/rlut_cfDay_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
tgcsw=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/rsutcs_cfDay_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
tgclw=ncassoc('/data/local2/hadyt/cmip5/'+models(m)+'/rlutcs_cfDay_'+models(m)+'_amip_r1i1p1_'+period1(m)+'-'+period2(m)+'.nc',/ignore)
for r=0,2 do begin
  cldarr=fltarr(49,nregimes(r))
  npts=fltarr(nregimes(r))
  ncf_clu=fltarr(nregimes(r))
  nncf_clu=fltarr(nregimes(r))

  snregimes=strcompress(string(nregimes(r)),/remove_all)
  reffile='auxiliary_data/global_2d_isccp_nclusters'+snregimes+'_'+refrgn(r)+'.data'

  openw,2,region(r)+'_onto_isccp.frac'
  openw,3,region(r)+'_onto_isccp.alb'
  openw,4,region(r)+'_onto_isccp.nctp'
  openw,5,region(r)+'_onto_isccp.clt'
  openw,6,region(r)+'_onto_isccp.swcf'
  openw,7,region(r)+'_onto_isccp.lwcf'
  openw,8,region(r)+'_onto_isccp.albf'
  openw,9,region(r)+'_onto_isccp.nswcf'
  openw,10,region(r)+'_onto_isccp.ncf'
  openw,11,region(r)+'_onto_isccp.nncf'
  openw,15,region(r)+'_onto_isccp.cltfrac'
  openw,16,region(r)+'_onto_isccp.swcffrac'
  openw,17,region(r)+'_onto_isccp.lwcffrac'
  openw,18,region(r)+'_onto_isccp.albffrac'
  openw,19,region(r)+'_onto_isccp.nswcffrac'
  openw,20,region(r)+'_onto_isccp.ncffrac'
  openw,21,region(r)+'_onto_isccp.nncffrac'

for mon=1,12 do begin
    print, 'mon',mon
    ; cloud regime calculation
    fuzzy_project_onto_clusters_2d_ncdf_mon, mon=mon, tgsic=tgsic, tgsnc=tgsnc, tgalbis=tgalbis, tgpctis=tgpctis, tgcltis=tgcltis, $
         reffile=reffile, refnseed=nregimes(r), regs=region(r), diagcode='cfmip', $
         /regrid, $
         tgrsdt=tgrsdt, tgnrsdt=tgnrsdt, tgtsw=tgtsw, tgtlw=tgtlw, tgcsw=tgcsw, tgclw=tgclw, $
         count_out=count_out, awo=awo, swcf_clu=swcf_clu, lwcf_clu=lwcf_clu, $
	 albf_clu=albf_clu, nswcf_clu=nswcf_clu, arr_out=arr_out
    ; net cloud radiative effect and annually normalized net cloud radiative effect
    ncf_clu=swcf_clu+lwcf_clu
    nncf_clu=nswcf_clu+lwcf_clu

    ; output results
    printf, format='(7(F10.3))',2,awo
    printf, format='(7(F10.3))',3,arr_out(0,0),arr_out(0,1),arr_out(0,2),arr_out(0,3),arr_out(0,4),arr_out(0,5),arr_out(0,6)
    printf, format='(7(F10.3))',4,arr_out(1,0),arr_out(1,1),arr_out(1,2),arr_out(1,3),arr_out(1,4),arr_out(1,5),arr_out(1,6)
    printf, format='(7(F10.3))',5,arr_out(2,0),arr_out(2,1),arr_out(2,2),arr_out(2,3),arr_out(2,4),arr_out(2,5),arr_out(2,6)
    printf, format='(7(F10.3))',6,swcf_clu
    printf, format='(7(F10.3))',7,lwcf_clu
    printf, format='(7(F10.3))',8,albf_clu
    printf, format='(7(F10.3))',9,nswcf_clu
    printf, format='(7(F10.3))',10,ncf_clu
    printf, format='(7(F10.3))',11,nncf_clu
    printf, format='(7(F10.3))',15,arr_out(2,*)*awo
    printf, format='(7(F10.3))',16,swcf_clu*awo
    printf, format='(7(F10.3))',17,lwcf_clu*awo
    printf, format='(7(F10.3))',18,albf_clu*awo
    printf, format='(7(F10.3))',19,nswcf_clu*awo
    printf, format='(7(F10.3))',20,ncf_clu*awo
    printf, format='(7(F10.3))',21,nncf_clu*awo
  
endfor
close,2
close,3
close,4
close,5
close,6
close,7
close,8
close,9
close,10
close,11
close,15
close,16
close,17
close,18
close,19
close,20
close,21
endfor
endfor
end
