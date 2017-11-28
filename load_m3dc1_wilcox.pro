; Load an M3D-C1 output file and dump the relevant data into an IDL .sav file
; This is not in a good state right now! Use something else.



save_fileloc = '/u/wilcoxr/m3dc1/170007_3625/test170007.sav'
h5_fileloc = 'equilibrium_l.h5'

; Do this first (in IDL interactive session):
;@colors_kc.pro
;.compile /u/shafer/idl/profdb/get_profdb.pro
;.compile /u/shafer/idl/addferraro.pro
;.compile /u/shafer/idl/quickpsinorm.pro
;addferraro

a=reada(157306,3750,runid='EFIT05')
g=readg('/u/shafer/efit/157306/g157306.03750_m3dc1')
psin2d = quickpsinorm(g)

g1 = readg(157306,3750)
n1z = 501
psin2d1 = quickpsinorm(g,nz=n1z,r=g1r,z=g1z)
n1r = n_elements(g1r)

;;restore,'/u/shafer/idlsav/besr_157306_80.sav'   ;besstr

dir='/u/ferraro/data/157306/03750/rw2_n=3/' & mtag = '2f'


getdat = 1
eqslice = -1
pslice = 1



file = dir+'C1.h5'
field = 'pe'
tag = '157306_3750_rw2_n=3'
shot = 157306
time = 3750
prunid='rec3d'

linfac = 2.
linfac = -2. ;why wasn't this set before!!!!
linfac = -2.5*4/!pi
ntor=3
pts = 401

nprof = 701
zv = fltarr(nprof)
rv = findgen(nprof)/(nprof-1)*0.3+2.1

rrange = [2,2.5]
zrange = [-0.5,0.5]
rrange = [1.9,2.5]
zrange = [-0.2,0.9]

fields = ['te','den','ti','pe','p']
labels = ['L!DTe!N','L!Dne!N','L!DTi!N','L!DPe!N','L!DP!N']
labels2= ['T!De!N','n!De!N','T!Di!N','P!De!N','P']

;fields = ['te','den']
;labels = ['L!DTe!N','L!Dne!N']
;labels2= ['T!De!N','n!De!N']


nfields = n_elements(fields)

phi_bes = 146.
phi_1 = (((30 - phi_bes) + 360) mod 360)
phi_2 = ((((30 - phi_bes) + 360) mod 360) + 60)

phis = [phi_1,phi_2]
nphis = n_elements(phis)

ntot = 360
dtot = 1 
phi_machine_tot = findgen(ntot)*dtot
;phi_m3dc1_tot = (((30. - phi_machine_tot) + 360) mod 360)
phi_m3dc1_tot = 30. - phi_machine_tot

phi_thom = 120

phit_1 = (((30 - phi_thom) + 360) mod 360)
phit_2 = ((((30 - phi_thom) + 360) mod 360) + 60)
phist = [phit_1,phit_2]
nphist = n_elements(phist)

;reflec
phi_ref = 255
phir_1 = (((30 - phi_ref) + 360) mod 360)
phir_2 = ((((30 - phi_ref) + 360) mod 360) + 60)
phisr = [phir_1,phir_2]
nphisr = n_elements(phisr)


nzt = 301
ztr = [-6.57,82.58]/1e2
zt = findgen(nzt)/(nzt-1)*(ztr[1]-ztr[0])+ztr[0]
rt = replicate(1.94,nzt)

if getdat eq 1 then begin

    p = get_profdb(shot,time,prunid)
    p1 = get_profdb(shot,time,prunid)
    p2 = get_profdb(shot,3850,'bw01')
    
    eqs = fltarr(nprof,nfields)
    perts = fltarr(nprof,nfields,nphis)
    perts_tot = fltarr(nprof,nfields,ntot)
    pertst = fltarr(nzt,nfields,nphist)
    pertsr = fltarr(nprof,nfields,nphisr)
    eqst = fltarr(nzt,nfields)


    for i = 0,nfields-1 do begin

        rese = read_field(fields[i],r,z,t,slice=eqslice,filename=file,points=pts,$
                          /mks,/linear,/equilibrium,rrange=rrange,zrange=zrange)

        int1 = interpol(lindgen(pts),r,rv)
        int2 = interpol(lindgen(pts),z,zv)
        eqs[*,i] = interpolate(reform(rese),int1,int2)

        int1 = interpol(lindgen(pts),r,rt)
        int2 = interpol(lindgen(pts),z,zt)
        eqst[*,i] = interpolate(reform(rese),int1,int2)

        res = read_field(fields[i],r,z,t,slice=pslice,filename=file,points=pts,/mks,$
                         /linear,/complex,rrange=rrange,zrange=zrange)

        imag = imaginary(reform(res))
        real = real_part(reform(res))   ;same as resreal...
        amp = linfac*sqrt(imag^2+real^2)
        phase = atan(reform(res),/phase)

        ;midplane
        int1 = interpol(lindgen(pts),r,rv)
        int2 = interpol(lindgen(pts),z,zv)
        
        for j=0,nphis-1 do begin
            
            perts[*,i,j] = interpolate(amp*cos(phase+phis[j]/!radeg*ntor),int1,int2)
        
        endfor


        for j=0,ntot-1 do begin
            
            perts_tot[*,i,j] = interpolate(amp*cos(phase+phi_m3dc1_tot[j]/$
                                            !radeg*ntor),int1,int2)
        
        endfor

        for j=0,nphisr-1 do begin
            
            pertsr[*,i,j] = interpolate(amp*cos(phase+phisr[j]/!radeg*ntor),int1,int2)
                    
        endfor

        ;thomson
        int1 = interpol(lindgen(pts),r,rt)
        int2 = interpol(lindgen(pts),z,zt)
        
        for j=0,nphist-1 do begin
            
            pertst[*,i,j] = interpolate(amp*cos(phase+phist[j]/!radeg*ntor),int1,int2)
                    
        endfor

        ;res_p1 = amp*cos(phase+phi_1/!radeg*ntor)
        ;res_p2 = amp*cos(phase+phi_2/!radeg*ntor)

        ;cut_eq = interpolate(rese,int1,int2)
        ;cut_p1 = interpolate(res_p1,int1,int2)
        ;cut_p2 = interpolate(res_p2,int1,int2)
        
    endfor

endif


;gradient scale lengths
eqsl = eqs*0
pertsl = perts*0
pertslr = pertsr*0
pertsl_tot = perts_tot*0


for i = 0,nfields-1 do begin
    
    eqsl[*,i] = a.d.aminor/1e2/(eqs[*,i]/deriv(rv,eqs[*,i]))
    
    for j=0,nphis-1 do begin
    
        val = eqs[*,i]+perts[*,i,j]
        ;pertsl[*,i,j] = a.d.aminor/1e2/(val/deriv(rv,val))
        ;part = 1./eqs[*,i]*deriv(rv,perts[*,i,j])-eqs[*,i]/perts[*,i,j]*deriv(rv,eqs[*,i])
        part = 1./eqs[*,i]*deriv(rv,perts[*,i,j])-perts[*,i,j]/eqs[*,i]^2*deriv(rv,eqs[*,i])
        pertsl[*,i,j] = a.d.aminor/1e2*part + eqsl[*,i]
    
    endfor

    for j=0,nphisr-1 do begin
    
        val = eqs[*,i]+pertsr[*,i,j]
        part = 1./eqs[*,i]*deriv(rv,pertsr[*,i,j])-pertsr[*,i,j]/eqs[*,i]^2*deriv(rv,eqs[*,i])
        pertslr[*,i,j] = a.d.aminor/1e2*part + eqsl[*,i]
    
    endfor


    for j=0,ntot-1 do begin
    
        val = eqs[*,i]+perts_tot[*,i,j]
        part = 1./eqs[*,i]*deriv(rv,perts_tot[*,i,j])-perts_tot[*,i,j]/$
               eqs[*,i]^2*deriv(rv,eqs[*,i])
        pertsl_tot[*,i,j] = a.d.aminor/1e2*part + eqsl[*,i]
    
    endfor


endfor


int1 = interpol(lindgen(g.mw),g.r,rv)
int2 = interpol(lindgen(g.mh),g.z,zv)
psinrv = interpolate(psin2d,int1,int2)



plotopt = 12


zr1 = [-1,1]*3.5e20
zr1= [-1,1]*3e18*1.6




    @greek_chars.pro
    @colors_kc.pro
    !p.background = white
    !p.color = black
    xr = [2.2,2.3]
    cols = [blue,red]
    cols = [red,blue]




if plotopt eq 12 then begin

    rad = psinrv & xr = [0.85,1.05] & sepval = 1 & xtitle = psinstr

    if !d.name eq 'PS' then !p.charsize=1.5
   ;------------------------------------
    ;   m3dc1 density perturbation in phi
    load_difct
    i = where(fields eq 'den') & pow = 3
    
    ;for i = 0,nfields-1 do begin
    
    val = reform(perts_tot[*,i,*])
    vald = val*0
    for j = 0, ntot-1 do begin
        vald[*,j] = deriv(rv,val[*,j]);+eqs[*,i])
    endfor
    zr1 = [-1,1]*3.5e20
    zr1= [-1,1]*3e18
    pixplotcb,val,psinrv,phi_machine_tot,xr=xr,/autoz,$
        ytitle=phistr+ ' (deg.)',xtitle=xtitle,cbarlabel=fields[i]
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    ;oplot,besstr.psf_avgs,besstr.phipos*!radeg,col=black
    ;oplot,besstr.psf_avgs,besstr.phipos*!radeg,col=black,psym=sym(1)

    oplot,[0,10],[1,1]*146,linestyle=2,thick=3
    oplot,[0,10],[1,1]*255,linestyle=2,thick=2

    ;endfor

    xyouts,0.2,0.01,file,/normal

    !p.multi=[0,4,1]
    !p.charsize = 1
endif



!p.charsize=1
!p.multi=[0,1,1]
!p.position=[0,0,0,0]
    
save,r,z,psinrv,eqs,perts,perts_tot,pertst,pertsr,eqst,pertsl_tot,fields,filename=save_fileloc


end
