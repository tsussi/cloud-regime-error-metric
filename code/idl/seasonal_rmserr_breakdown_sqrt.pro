;Estimating breakdown of errors in [Wm-2] unit

regs=['2020','exts','extn']
nregs=n_elements(regs)
;for r=0,nregs-1 do begin
for r=0,0 do begin
if (r eq 0) then begin
labels=['Shal. Cu.','Congest.','Thin Cirr.','Trans.','Anv. Cirr.','Deep Conv.','Stratocu.']
endif
if (r gt 0) then begin
labels=['Shal. Cu.','Congest.','Trans.','Cirrus','Stratocu.','Frontal','Thin Cirr.']
endif

expt='mon'
obs=['auxiliary_data']
models=['auxiliary_data','cmip5/HadGEM2-A','cmip5/CanAM4','cmip5/MRI-CGCM3', 'cmip5/MIROC5','cmip5/IPSL-CM5B-LR','cmip5/CNRM-CM5']

obs_names=['ISCCP']
vars=['frac','swcf','nswcf','ncf','nncf','lwcf','clt','swcffrac','nswcffrac','ncffrac','nncffrac','lwcffrac','cltfrac']

nmonths=12
nmodels=n_elements(models)
nregimes=n_elements(labels)
nvars=n_elements(vars)

;output climatological variables
for v=0,nvars-1 do begin
  openw,10,'rmserr_breakdown_sqrt_'+vars(v)+'_'+regs(r)+'_onto_isccp'
for m=2,2 do begin
   openr,1,obs(0)+'/'+expt+'/'+regs(r)+'_onto_isccp.'+vars(v)+''
   openr,11, regs(r)+'_onto_isccp.'+vars(v)+''
;count -> frac conversion which is made beforehand
  vo=fltarr(nregimes,nmonths)
  vm=fltarr(nregimes,nmonths)
  climo=fltarr(nregimes)
  climm=fltarr(nregimes)
  varo=fltarr(nregimes)
  varm=fltarr(nregimes)
  varclimmo=fltarr(nregimes)
  dstdev=fltarr(nregimes)
  dstdev2=fltarr(nregimes)
  varanommo=fltarr(nregimes)
  covterm=fltarr(nregimes)
  varphase=fltarr(nregimes)
  covphase=fltarr(nregimes)
  varms=fltarr(nregimes,nmodels)
;frac
for mon=0,nmonths-1 do begin
  tmp1=fltarr(nregimes)
  readf,1,tmp1
  vo(*,mon)=tmp1
  readf,11,tmp1
  vm(*,mon)=tmp1
endfor
close,1
close,11
;climatological annual mean
  climo(*)=(vo(*,0)+vo(*,1)+vo(*,2)+vo(*,3)+vo(*,4)+vo(*,5)+vo(*,6)+vo(*,7)+vo(*,8)+vo(*,9)+vo(*,10)+vo(*,11))/12.
  climm(*)=(vm(*,0)+vm(*,1)+vm(*,2)+vm(*,3)+vm(*,4)+vm(*,5)+vm(*,6)+vm(*,7)+vm(*,8)+vm(*,9)+vm(*,10)+vm(*,11))/12.
print, 'climo(*)',climo(*)
print, 'climm(*)',climm(*)
  varclimmo(*)=(climm(*)-climo(*))^2
for mon=0,nmonths-1 do begin
;variance
;Ao**2
  varo(*)=varo(*)+(vo(*,mon)-climo(*))^2/12.
;Am**2
  varm(*)=varm(*)+(vm(*,mon)-climm(*))^2/12.
;(Am-Ao)**2
  dstdev(*)=sqrt(varm(*))-sqrt(varo(*))
  dstdev2(*)=(sqrt(varm(*))-sqrt(varo(*)))^2
;Centered rmserr: Variance of the difference in model and observation's amonthly
;anomaly from each climatology.
  varanommo(*)=varanommo(*)+((vm(*,mon)-climm(*))-(vo(*,mon)-climo(*)))^2/12.
endfor
;Obtain covariance term by subtracting (Am-Ao)^2)
  covterm(*)=varanommo(*)-dstdev2(*)
  for i=0,nregimes-1 do begin
  if (covterm(i) lt 0) then covterm(i)=0.
  endfor
  varphase(*)=covterm(*)*sqrt(varo(*))/sqrt(varm(*))
  covphase(*)=covterm(*)-varphase(*)
;sqrt
  varclimmo(*)=sqrt(varclimmo(*))
  varo(*)=sqrt(varo(*))
  varm(*)=sqrt(varm(*))
  varanommo(*)=sqrt(varanommo(*))
  covterm(*)=sqrt(covterm(*))
  varphase(*)=sqrt(varphase(*))
  covphase(*)=sqrt(covphase(*))

;write to a file
if (m gt 0) then printf, format='(14(F10.4))',10,dstdev(0),covterm(0),dstdev(1),covterm(1),dstdev(2),covterm(2),dstdev(3),covterm(3),dstdev(4),covterm(4),dstdev(5),covterm(5),dstdev(6),covterm(6)
endfor
close,10
endfor
endfor
end
