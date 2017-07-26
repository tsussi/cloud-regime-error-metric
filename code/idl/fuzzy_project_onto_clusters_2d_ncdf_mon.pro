; NAME:fuzzy_project_onto_clusters_2d_ncdf_mon
; PURPOSE: Project ISCCP diagnostics onto a set of cluster histograms using the
;          lightweight ISCCP mean CTP, tau and TTC.
; REQUIRED INPUTS:
;    mon - Month of data which are processed
;    tgsic, tgsnc - Mask for ice/snow covered area
;    tgalbis, tgpctis, tgcltis - ISCCP diagnostics
;    reffile - file containing reference clusters to be projected onto
;    refnseed - number of clusters in reference file
;    regs - region: can be one of '2020','exts','extn',icey','other' if other then region_area must be set
;    diagcode - identifies how isccp diagnostics are ordered. Can be one of: 
;               'isccp' (use for data in /project/earthobs)
;               'cfmip' (use for data in /project/cfmip)
;               'um' (use for data from the UM) 
; OPTIONAL INPUTS:
;    region_area - co-ordinates of area to be extracted for non-standard regions (required if regs='other')
;    tgrsdt, tgnrsdt - TOA incoming SW flux (required if tge is set). If not equal to tge then an annual cycle is required
;    tgtsw, tgtlw, tgcsw, tgclw - TOA flux diagnostics (if unset then assumed unavailable)
; KEYWORD PARAMETERS:
;    /isccp_ice - Daily ISCCP-IS data being used for the ice mask
;    /regrid - regrid onto ISCCP 2.5 degree grid
;    /erbe - clear-sky fluxes not avail when grid-point is cloud covered so use monthly means
; OUTPUTS
;    count_out(refnseed) - number of datapoints assigned to each cluster
;    awo(refnseed) -      area-weighted occurence of each cluster
;    swcf_clu(refnseed) - average shortwave cloud forcing for each cluster
;    lwcf_clu(refnseed) - average longwave cloud forcing for each cluster
;    nswcf_clu(refnseed) - average SWCF normalised by the insolation 
;    arr_out(3,refnseed) - average albedo, CTP, TCC for each cluster
; NOTES:
;    Assumes that none of the reference clusters are clear-sky

pro fuzzy_project_onto_clusters_2d_ncdf_mon, mon=m, tgsic=tgsic, tgsnc=tgsnc, tgalbis=tgalbis, tgpctis=tgpctis, tgcltis=tgcltis, $
         reffile=reffile, refnseed=refnseed, regs=regs, diagcode=diagcode, $
	 isccp_ice=isccp_ice, regrid=regrid, erbe=erbe, region_area=region_area, $
         tgrsdt=tgrsdt, tgnrsdt=tgnrsdt, tgtsw=tgtsw, tgtlw=tgtlw, tgcsw=tgcsw, tgclw=tgclw, $
         count_out=count_out, awo=awo, swcf_clu=swcf_clu, lwcf_clu=lwcf_clu, $
	 albf_clu=albf_clu, nswcf_clu=nswcf_clu, arr_out=arr_out
	 
@math_startup
@stat_startup
crf_avail=1
if not(keyword_set(tgsnc)) then ice_req=0 else ice_req=1
if (regs ne 'other') then region_area=[0,-90,360,90]
nseeds=strcompress(string(refnseed),/remove_all)

dummy=fltarr(refnseed)
ref_work=fltarr(3,refnseed)
ref_groups=fltarr(3,refnseed)

; Reference groups
openr,1,reffile
readf,1,dummy
readf,1,ref_work
readf,1,dummy
readf,1,dummy
readf,1,dummy
close,1

ref_groups=ref_work

if keyword_set(regrid) then begin
  tgregrid=ncassoc('auxiliary_data/cltisccp_obs4MIPs_ISCCP_L3_V1.0_200511-200511.nc',/ignore)
  get_size_glob=pp_extract(ppa(tgregrid,ss(from=[2005,11],to=[2005,11]),/all),region_area)
endif

get_ndays=ppa(tgcltis,ss(from=[1979,1,1],to=[1979,12,31],monv=m),/all)
ndays=n_elements(get_ndays)
s=size(get_size_glob.data)

print, 'number of data times: ', ndays

lats = pp_coords(get_size_glob, /ygen)

;Create a mask for ice/snow covered area
sic=ppa(tgsic, ss(from=[1979,1,1],to=[1979,12,31]),/all)
snc=pp_regrid(ppa(tgsnc, ss(from=[1979,1,1],to=[1979,12,31]),/all),sic)
sicmask=pp_avg(sic)
sncmask=pp_avg(snc)
sicsncmask=sicmask
imdiis=where(sicmask.data gt 0. or sncmask.data gt 0.)
sicsncmask.data(imdiis)=1.
imdinois1=where(sicmask.data eq 0. and sncmask.data eq 0.)
sicsncmask.data(imdinois1)=sicmask.bmdi
imdinois2=where(sicmask.data eq 0. and sncmask.data eq sncmask.bmdi,count)
if (count gt 0) then sicsncmask.data(imdinois2)=sicmask.bmdi
tgice=sicsncmask

if (ice_req eq 1) then begin
  if (diagcode eq 'isccp' or keyword_set(isccp_ice)) then begin
    isice=pp_regrid(ppa(tgice, ss(user4=99931),/all),get_size_glob)
    if (keyword_set(isccp_ice)) then begin
      ice=replicate(isice(0),ndays)
      for i=0,ndays-1 do begin
        pt=where(isice.lbmon eq get_ndays(i).lbmon AND isice.lbdat eq get_ndays(i).lbdat, count)
        if (count eq 0) then begin
          pt=where(isice.lbmon eq get_ndays(i).lbmon and isice.lbdat eq get_ndays(i).lbdat-1, count)
          if (count eq 0) then begin
            pt=where(isice.lbmon eq get_ndays(i).lbmon and isice.lbdat eq get_ndays(i).lbdat-2, count)
          endif
        endif
        ice(i)=isice(pt(0))
      endfor
    endif else begin
      ice=isice
    endelse
  endif else begin
    isice=pp_regrid(ppa(tgice, /all),get_size_glob)
    ice=replicate(isice(0),ndays)
    for i=0,ndays-1 do begin
      pts=where(ice(i).data eq ice(i).bmdi)
      ice(i).data(pts)=-1.0E30
      ice(i).bmdi=-1.0E30
    endfor
  endelse
  mask=ice
  mask.data=mask(0).bmdi

  tarea=0.
  countall=0
  for i=0,ndays-1 do begin
     if (regs eq '2020') then begin
        pts=where(lats ge -20 and lats le 20,count)
        mask(i).data(*,pts)=1.0
        tarea=tarea+pp_area_total(mask)
     endif 
     if (regs eq 'exts') then begin
        pts=where(lats lt -20)
        mask(i).data(*,pts)=1.0
;simask
        pts=where(ice(i).data gt 0.0)
        mask(i).data(pts)=mask(0).bmdi
        pts=where(mask(i).data eq 1.0,count)
        tarea=tarea+pp_area_total(mask)
if (i eq 0) then print, 'regs,after', regs,mask(0).data
     endif
     if (regs eq 'extn') then begin
        pts=where(lats gt 20)
        mask(i).data(*,pts)=1.0
;simask
        pts=where(ice(i).data gt 0.0)
        mask(i).data(pts)=mask(0).bmdi
        pts=where(mask(i).data eq 1.0,count)
        tarea=tarea+pp_area_total(mask)
     endif
     if (regs eq 'icey') then begin
        pts=where(ice(i).data le 0.0,count)
        mask(i).data(pts)=mask(0).bmdi
        tarea=tarea+pp_area_total(mask)
     endif
    countall=countall+count
  endfor
endif else begin
  mask=replicate(get_size_glob,ndays)
  mask.data=1.0
endelse
print, 'Countall',countall
index=where(mask(*).data eq 1.0, count)
print, 'Count',count
tauctp_data=fltarr(3,count)
avg_mask=pp_avg(mask,mdtol=1)

;Process isccp diagnostics
print, 'Reading ISCCP diagnostics'

if (diagcode eq 'isccp' or diagcode eq 'cfmip') then begin
  ppfield_glob=pp_extract(pp_regrid(ppa(tgcltis,ss(from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob),region_area)
  ppfield_clt=ppfield_glob
  ppfield_glob=pp_ff('a/100.',ppfield_glob,/quiet)
  ppfield_data=ppfield_glob(*).data
  tauctp_data(2,*)=ppfield_data(index)
  ppfield_glob=pp_extract(pp_regrid(ppa(tgalbis,ss(from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob),region_area)
  ppfield_alb=ppfield_glob
  ppfield_glob=pp_ff('a/b*100.',ppfield_alb, ppfield_clt,/quiet)
  ppfield_data=ppfield_glob(*).data
  tauctp_data(0,*)=ppfield_data(index)
  ppfield_glob=pp_extract(pp_regrid(ppa(tgpctis,ss(from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob),region_area)
  ppfield_pct=ppfield_glob
  ppfield_glob=pp_ff('a/b*100.',ppfield_pct, ppfield_clt,/quiet)
;cloud top pressure in Pa.
  ppfield_glob=pp_ff('a/100000.',ppfield_glob,/quiet)
  ppfield_data=ppfield_glob(*).data
  tauctp_data(1,*)=ppfield_data(index)
endif else begin
  weight=ppa(tgi, ss(st=2330),/all)
  cld_weight=ppa(tgi, ss(st=2334),/all)
  ppfield_glob=ppa(tgi, ss(user4=2331),/all)
  ppfield_glob=pp_regrid(pp_ff('a/b',ppfield_glob,cld_weight,/math_fix,/quiet),get_size_glob)
  ppfield_data=ppfield_glob(*).data
  tauctp_data(0,*)=ppfield_data(index)
  ppfield_glob=ppa(tgi, ss(user4=2333),/all)
;cloud top pressure in Pa.
  ppfield_glob=pp_regrid(pp_ff('a/(100000.*b)',ppfield_glob,cld_weight,/math_fix,/quiet),get_size_glob)
  ppfield_data=ppfield_glob(*).data
  tauctp_data(1,*)=ppfield_data(index)
  ppfield_glob=pp_regrid(pp_ff('a/b',cld_weight,weight,/math_fix,/quiet),get_size_glob)
  ppfield_data=ppfield_glob(*).data
  tauctp_data(2,*)=ppfield_data(index)
endelse
taumdi=ppfield_glob(0).bmdi

if (crf_avail eq 1) then begin
  outsw=pp_extract(pp_regrid(ppa(tgtsw,ss(from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob),region_area)
  outlw=pp_extract(pp_regrid(ppa(tgtlw,ss(from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob),region_area)
  if keyword_set(erbe) then begin
    start=intarr(3)
    finish=intarr(3)
    start(0)=outsw(0).lbyr
    start(1)=outsw(0).lbmon
    start(2)=outsw(0).lbdat
    finish(0)=outsw(ndays-1).lbyr
    finish(1)=outsw(ndays-1).lbmon
    finish(2)=outsw(ndays-1).lbdat
    csoutsw=outsw
    csuplw=outlw
    for y=start(0),finish(0) do begin
      if (y eq start(0)) then mstart=start(1) else mstart=1
      if (y eq finish(0)) then mfinish=finish(1) else mfinish=12
      for m=mstart,mfinish do begin
        swpts=where((outsw.lbyr eq y) AND (outsw.lbmon eq m))
        lwpts=where((outlw.lbyr eq y) AND (outlw.lbmon eq m))
        cssw=pp_regrid(ppa(tgcsw,ss(st=1209,from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob)
        cslw=pp_regrid(ppa(tgclw,ss(st=2206,from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob)
        csoutsw(swpts)=pp_avg(cssw,mdtol=1)
        csuplw(swpts)=pp_avg(cslw,mdtol=1)
      endfor
    endfor
  endif else begin
    csoutsw=pp_extract(pp_regrid(ppa(tgcsw,ss(from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob),region_area)
    csuplw=pp_extract(pp_regrid(ppa(tgclw,ss(from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob),region_area)
  endelse
  swcf=pp_diff(csoutsw,outsw)
  lwcf=pp_diff(csuplw,outlw)

;Assume all models submit rsdt  
  solar=pp_extract(pp_regrid(ppa(tgrsdt,ss(from=[1979,1,1],to=[1979,12,31],monv=m),/all),get_size_glob),region_area)
  nsolarin=pp_avg(pp_regrid(ppa(tgnrsdt,ss(from=[1979,1,1],to=[1979,12,31]),/all),get_size_glob))
  solaravg=replicate(nsolarin,ndays)
  albf=pp_ff('-a/b',swcf,solar, /math_fix)
  nswcf=pp_ff('-a*b',albf,solaravg)

  swcfmdi=swcf(0).bmdi
  lwcfmdi=lwcf(0).bmdi

  swcf_data=swcf.data
  lwcf_data=lwcf.data
  albf_data=albf.data
  nswcf_data=nswcf.data
  solar_data=solar.data
  swcf_data=swcf_data(index)
  lwcf_data=lwcf_data(index)
  albf_data=albf_data(index)
  nswcf_data=nswcf_data(index)
  pts=where(nswcf_data ne swcfmdi)
endif

area=replicate(get_size_glob, ndays)
areamdi=area(0).bmdi
for i=0,ndays-1 do begin 
  if keyword_set(area_not_avail) then begin
    area(i).data=1.0
  endif else begin
    area(i).data=pp_area(area(i))
  endelse
endfor

area_data=area.data
area_data=area_data(index)

print, 'Projecting data points'

points=where((tauctp_data(0,*) NE taumdi) AND (tauctp_data(1,*) NE taumdi) AND (tauctp_data(2,*) NE taumdi) AND (area_data NE areamdi), fred)

group=fltarr(fred)
count_out=fltarr(refnseed)
awo=fltarr(refnseed)
area_r=fltarr(refnseed)
swcf_clu=fltarr(refnseed)
lwcf_clu=fltarr(refnseed)
albf_clu=fltarr(refnseed)
nswcf_clu=fltarr(refnseed)
arr_out=fltarr(3,refnseed)
locate=replicate(get_size_glob,refnseed)

location_map=replicate(get_size_glob,ndays)

for i=0L,fred-1 do group(i)=min_euclidean(tauctp_data(*,points(i)),ref_groups(*,0:refnseed-1))

for i=0,refnseed-1 do begin
  mem=where(group eq (i), count) 
  if (crf_avail eq 1) then begin
    swcfmem=where(group eq (i) AND swcf_data(points) ne swcfmdi) 
    lwcfmem=where(group eq (i) AND lwcf_data(points) ne lwcfmdi) 
    albfmem=where(group eq (i) AND albf_data(points) ne swcfmdi)
    nswcfmem=where(group eq (i) AND nswcf_data(points) ne swcfmdi)
  endif 
  if (count gt 0) then begin
    count_out(i)=count
    for j=0,2 do begin
       arr_out(j,i)=avg(tauctp_data(j,points(mem))*area_data(points(mem)))/avg(area_data(points(mem)))
    endfor
    if (crf_avail eq 1) then begin
      swcf_clu(i)=avg(swcf_data(points(swcfmem))*area_data(points(swcfmem)))/avg(area_data(points(swcfmem)))
      lwcf_clu(i)=avg(lwcf_data(points(lwcfmem))*area_data(points(lwcfmem)))/avg(area_data(points(lwcfmem)))
      albf_clu(i)=avg(albf_data(points(albfmem))*area_data(points(albfmem)))/avg(area_data(points(albfmem)))
      nswcf_clu(i)=avg(nswcf_data(points(nswcfmem))*area_data(points(nswcfmem)))/avg(area_data(points(nswcfmem)))
      awo(i)=avg(area_data(points(mem)))*count/tarea(0)
    endif 
  endif else begin
    count_out(i)=0
    arr_out(*,i)=999.99  
    swcf_clu(i)=999.99
    lwcf_clu(i)=999.99
    albf_clu(i)=999.99
    nswcf_clu(i)=999.99
  endelse
  count_out(i)=count_out(i)/countall

endfor
locate(*).data=0.0
data_cut=0
fdata_start=0
for i=0,ndays-1 do begin
  workmap=get_size_glob
  workmap.lbyr=get_ndays(i).lbyr
  workmap.lbmon=get_ndays(i).lbmon
  workmap.lbdat=get_ndays(i).lbdat
  workmap.lbday=get_ndays(i).lbday
  workmap.lbyrd=get_ndays(i).lbyrd
  workmap.lbmond=get_ndays(i).lbmond
  workmap.lbdatd=get_ndays(i).lbdatd
  workmap.lbdayd=get_ndays(i).lbdayd
  data_start=i*s(1)*s(2)
  data_end=((i+1)*s(1)*s(2))-1
  workmap.data=workmap.bmdi
  pts=where(index ge data_start and index le data_end, count)
  fdata_end=fdata_start+count
  dindex=index(pts)-data_start
  refill=where((tauctp_data(0,fdata_start:fdata_end-1) NE taumdi) AND(tauctp_data(1,fdata_start:fdata_end-1) NE taumdi) AND(tauctp_data(2,fdata_start:fdata_end-1) NE taumdi) AND (area_data(fdata_start:fdata_end-1) NE areamdi),fred)
  if (fred gt 0) then begin
    workmap.data(dindex(refill))=group(data_cut:data_cut+fred-1)
    location_map(i)=pp_ff('a+2',workmap)
    for j=0,refnseed-1 do begin
      pts=where(workmap.data eq j, count)
      if (count gt 0) then locate(j).data(pts)=locate(j).data(pts)+1.0/ndays
    endfor
  endif
  data_cut=data_cut+fred
  fdata_start=fdata_end
endfor

pts=where(avg_mask.data eq avg_mask.bmdi, count)
if (count gt 0) then for i=0,refnseed-1 do locate(i).data(pts)=locate(0).bmdi
print, 'here_end'
if (count gt 0) then for i=0,refnseed-1 do count_out(i)=count_out(i)
return
end
