; adds a two column option
;
;addferraro

a=reada(157306,3750,runid='EFIT05')
g=readg('/u/shafer/efit/157306/g157306.03750_m3dc1')
psin2d = quickpsinorm(g)

g1 = readg(157306,3750)
n1z = 501
psin2d1 = quickpsinorm(g,nz=n1z,r=g1r,z=g1z)
n1r = n_elements(g1r)

restore,'/u/shafer/idlsav/besr_157306_80.sav'   ;besstr

dir='/u/ferraro/data/157306/03750/rw2_n=3/' & mtag = '2f'
;dir='/u/ferraro/data/157306/03750/rw2_n=3_1f/' & mtag = '1f'
;dir = '/u/ferraro/data/157306/03750/rw2_n=3_1f_vC/'& mtag = '1f_vc'
;dir = '/u/ferraro/data/157306/03750/rw3_n=3_vC/' & mtag = '2f_vc'
;dir = '/u/ferraro/data/157306/03750/rw2_n=3_vX/' & mtag = '2f_vx'
;dir = '/u/ferraro/data/157306/03750/rw2_n=3_2p/' & mtag = '2f_2p'

;many directory loop.. to setup 
flutags = ['1f','2f']
vtags = ['slow','kinda_slow','measured','kinda_fast','fast'];+'_ExB'
bdir = '/u/ferraro/data/157306/03750/rw2_n=3_'
iii=0
jjj=0

;for iii=0,n_elements(flutags)-1 do begin
;for jjj=0,n_elements(vtags)-1 do begin

;mtag = flutags[iii]+'_'+vtags[jjj]
;dir = bdir+mtag+'/'

getdat = 1*0
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

fields = ['te','den','ti','p']
labels = ['L!DTe!N','L!Dne!N','L!DTi!N','L!DP!N']
labels2= ['T!De!N','n!De!N','T!Di!N','P!D!N']


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

    p=get_profdb(shot,time,prunid)
    p1=get_profdb(shot,time,prunid)
    p2=get_profdb(shot,3850,'bw01')

    ;reflectometer...
    mdsconnect,'atlas.gat.com'
    mdsopen,'electrons',shot
    refd=mdsvalue('\reflect_3dr')
    reft=mdsvalue('dim_of(\reflect_3dr,0)')
    refr=mdsvalue('dim_of(\reflect_3dr,1)')
    mdsclose
    mdsdisconnect

    avtim = 35
    times = [3760,3860]
    profiles = bag_thomson_get_profiles(shot,times,avtim,/elm)
    stats = BAG_THOMSON_GET_STATS(profiles,ENSEMBLE=ensemble)  
    tags = TAG_NAMES(profiles)
    ntags = n_elements(tags)

 
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

;entropy.... P/n^(5/3) .... not working
enl = fltarr(nprof,nphis)
ipe = where(fields eq 'pe')
iden = where(fields eq 'den')
ee = (1e19)^(3./5)*reform(eqs[*,ipe]/(eqs[*,iden]/1e19)^(5./3))
;ee = reform(eqs[*,ipe]/(eqs[*,iden])^(5./3))
enle = a.d.aminor/1e2/(ee/deriv(rv,ee))
for j=0,nphis-1 do begin
    
    pp = (1e19)^(3./5)*reform(perts[*,ipe,j]/(perts[*,iden,j]/1e19)^(5./3))
    ;pp = reform(perts[*,ipe,j]/(perts[*,iden,j])^(5./3))
    part = 1./ee*deriv(rv,pp)-pp/ee^2*deriv(rv,ee)
    ;enl[*,j]= a.d.aminor/1e2/(ee*deriv(rv,pp)-pp/ee^2*deriv(rv,ee))+enle
    enl[*,j] = a.d.aminor/1e2*part + enle
endfor


int1 = interpol(lindgen(g.mw),g.r,rv)
int2 = interpol(lindgen(g.mh),g.z,zv)
psinrv = interpolate(psin2d,int1,int2)

;int1 = interpol(lindgen(g.mh),g.z,0)
;rv_sep = interpol(g.bdry[0,*],g.bdry[1,*],0)

rv_sep = interpol(rv,psinrv,1)

int1 = interpol(lindgen(g.mw),g.r,rt)
int2 = interpol(lindgen(g.mh),g.z,zt)
psinzt = interpolate(psin2d,int1,int2)

zt_sep = interpol(zt,psinzt,1)


cinfo='host=huez.gat.com port=5681 user=d3dpub password=d3dpub dbname=d3d'
qres = pgsql_query("select timeid,runid,adjust_zts FROM profile_runs WHERE shot = 157306",connect_info=cinfo)
gqr = [where(qres.timeid eq 3750 and qres.runid eq 'rec3d'),where(qres.timeid eq 3850)]


nmap =n_elements(p1.map.r)
;int1 = interpol(lindgen(nzt),zt,p1.map.z);,/spline)
;int2 = interpol(lindgen(nzt),rt,p1.map.r);,/spline)
int1 = interpol(lindgen(nmap),p1.map.r,rt)
int2 = interpol(lindgen(nmap),p1.map.z,zt)


int1 = interpol(lindgen(g.mw),g.r,rt)
int2 = interpol(lindgen(g.mh),g.z,zt)


;p1psin = interpolate(p1.map.psirzn,int1,int2);,cubic=1)
p1psin = interpolate(psin2d,int1,int2);,cubic=1)

exshift = 0.01+.001;*0

p1te = interpol(p1.psi.tanh.t_e,p1.psi.tanh.psi_te,p1psin+qres[gqr[0]].adjust_zts-exshift)
p1ne = 10*interpol(p1.psi.tanh.n_e,p1.psi.tanh.psi_ne,p1psin+qres[gqr[0]].adjust_zts-exshift)

;nmap =n_elements(p2.map.r)
;int1 = interpol(lindgen(nmap),p2.map.r,rt)
;int2 = interpol(lindgen(nmap),p2.map.z,zt)
;p2psin = interpolate(p1.map.psirzn,int1,int2);,cubic=1)
p2psin = p1psin

p2te = interpol(p2.psi.tanh.t_e,p2.psi.tanh.psi_te,p2psin+qres[gqr[1]].adjust_zts-exshift)
p2ne = 10*interpol(p2.psi.tanh.n_e,p2.psi.tanh.psi_ne,p2psin+qres[gqr[1]].adjust_zts-exshift)

;nt1 = n_elements(p1.psi.dat.t_e)
;nt2 = n_elements(p2.psi.dat.t_e)

;p1th_z = fltarr(nt1)
;p2th_z = fltarr(nt2)

;p1th_z = interpol(zt,psinzt,p1.psi.dat.psi_te+qres[gqr[0]].adjust_zts-exshift)
;p2th_z = interpol(zt,psinzt,p2.psi.dat.psi_te+qres[gqr[1]].adjust_zts-exshift)

p1th_z = interpol(zt,psinzt,p1.psi.dat.psi_te-qres[gqr[0]].adjust_zts+exshift)
p2th_z = interpol(zt,psinzt,p2.psi.dat.psi_te-qres[gqr[1]].adjust_zts+exshift)

hbinsize = 0.0026 ; in meters
h1 = histogram(p1th_z,binsize=hbinsize,locations=loc)
nloc = n_elements(loc)

p1thzh = fltarr(nloc)
p1thzhe = fltarr(nloc)
p1teh = fltarr(nloc)
p1tehe = fltarr(nloc)
p1neh = fltarr(nloc)
p1nehe = fltarr(nloc)
p1thph = fltarr(nloc)

for i=0,nloc-1 do begin

    ind = where(p1th_z gt loc[i]-hbinsize/2. and p1th_z le loc[i]+hbinsize/2.)
    p1thzh[i] = mean(p1th_z[ind])
    p1thzhe[i] = stddev(p1th_z[ind])
    p1teh[i] = mean(p1.psi.dat.t_e[ind])
    p1tehe[i] = stddev(p1.psi.dat.t_e[ind])
    p1neh[i] = mean(p1.psi.dat.n_e[ind])*10
    p1nehe[i] = stddev(p1.psi.dat.n_e[ind])*10
    p1thph[i] = interpolate(psin2d,interpol(lindgen(g.mw),g.r,1.94),interpol(lindgen(g.mh),g.z,p1thzh[i]))

endfor

h2 = histogram(p2th_z,binsize=hbinsize,locations=loc)
nloc = n_elements(loc)

p2thzh = fltarr(nloc)
p2thzhe = fltarr(nloc)
p2teh = fltarr(nloc)
p2tehe = fltarr(nloc)
p2neh = fltarr(nloc)
p2nehe = fltarr(nloc)
p2thph = fltarr(nloc)

for i=0,nloc-1 do begin

    ind = where(p2th_z gt loc[i]-hbinsize/2. and p2th_z le loc[i]+hbinsize/2.)
    p2thzh[i] = mean(p2th_z[ind])
    p2thzhe[i] = stddev(p2th_z[ind])
    p2teh[i] = mean(p2.psi.dat.t_e[ind])
    p2tehe[i] = stddev(p2.psi.dat.t_e[ind])
    p2neh[i] = mean(p2.psi.dat.n_e[ind])*10
    p2nehe[i] = stddev(p2.psi.dat.n_e[ind])*10
    p2thph[i] = interpolate(psin2d,interpol(lindgen(g.mw),g.r,1.94),interpol(lindgen(g.mh),g.z,p2thzh[i]))


endfor

ravtim = avtim
;ravtim = 10

;reflectometer
win1 = where(reft gt times[0]-ravtim and reft lt times[0]+ravtim)
win2 = where(reft gt times[1]-ravtim and reft lt times[1]+ravtim)
hbinsize = 0.002 ; in meters
h1 = histogram(refr[win1,*],binsize=hbinsize,locations=loc)
nloc = n_elements(loc)

pr1r = fltarr(nloc)
pr1d = fltarr(nloc)
pr1re = fltarr(nloc)
pr1de = fltarr(nloc)
pr1psin = fltarr(nloc)
pr1psine = fltarr(nloc)


for i=0,nloc-1 do begin

    ind = where(refr[win1,*] gt loc[i]-hbinsize/2. and refr[win1,*] le loc[i]+hbinsize/2.)
    pr1r[i] = mean((refr[win1,*])[ind])
    pr1re[i] = stddev((refr[win1,*])[ind])
    pr1d[i] = mean((refd[win1,*])[ind])
    pr1de[i] = stddev((refd[win1,*])[ind])
    int1 = interpol(lindgen(g.mw),g.r,pr1r[i]+[0,pr1re[i],pr1re[i]*(-1)])
    int2 = interpol(lindgen(g.mh),g.z,replicate(0,3))
    tmp = interpolate(psin2d,int1,int2)
    ;int1 = interpol(lindgen(n1r),g1r,pr1r[i]+[0,pr1re[i],pr1re[i]*(-1)])
    ;int2 = interpol(lindgen(n1z),g1z,replicate(0,3))
    ;tmp = interpolate(psin2d1,int1,int2)

    pr1psin[i] = tmp[0]
    pr1psine[i] = abs(tmp[2]-tmp[1])
    
endfor

h2 = histogram(refr[win2,*],binsize=hbinsize,locations=loc)
nloc = n_elements(loc)

pr2r = fltarr(nloc)
pr2d = fltarr(nloc)
pr2re = fltarr(nloc)
pr2de = fltarr(nloc)
pr2psin = fltarr(nloc)
pr2psine = fltarr(nloc)

for i=0,nloc-1 do begin

    ind = where(refr[win2,*] gt loc[i]-hbinsize/2. and refr[win2,*] le loc[i]+hbinsize/2.)
    pr2r[i] = mean((refr[win2,*])[ind])
    pr2re[i] = stddev((refr[win2,*])[ind])
    pr2d[i] = mean((refd[win2,*])[ind])
    pr2de[i] = stddev((refd[win2,*])[ind])
    int1 = interpol(lindgen(g.mw),g.r,pr2r[i]+[0,pr2re[i],pr2re[i]*(-1)])
    int2 = interpol(lindgen(g.mh),g.z,replicate(0,3))
    tmp = interpolate(psin2d,int1,int2)
    ;int1 = interpol(lindgen(n1r),g1r,pr2r[i]+[0,pr2re[i],pr2re[i]*(-1)])
    ;int2 = interpol(lindgen(n1z),g1z,replicate(0,3))
    ;tmp = interpolate(psin2d1,int1,int2)

    pr2psin[i] = tmp[0]
    pr2psine[i] = abs(tmp[2]-tmp[1])
    
endfor





;for i = 0,nt1-1 do begin
;    
;    p1th_z[i] = interpol(zt,psinzt,p1.psi.dat.psi_te[i]+qres[gqr[0]].adjust_zts-exshift)
;        
;endfor

;for i = 0,nt2-1 do begin
;    
;    p2th_z[i] = interpol(p2.psi.dat.psi_te[i]+qres[gqr[1]].adjust_zts-exshift,psinzt,zt)
;        
;endfor



;p1ted = interpol(p1.psi.dat.t_e,p1.psi.dat.psi_te,p1psin+qres[gqr[0]].adjust_zts-exshift)
;p1ned = 10*interpol(p1.psi.dat.n_e,p1.psi.dat.psi_ne,p1psin+qres[gqr[0]].adjust_zts-exshift)

;p2ted = interpol(p2.psi.dat.t_e,p2.psi.dat.psi_te,p2psin+qres[gqr[1]].adjust_zts-exshift)
;p2ned = 10*interpol(p2.psi.dat.n_e,p2.psi.dat.psi_ne,p2psin+qres[gqr[1]].adjust_zts-exshift)



    dp_shift = 0.001 ;bob wants to see what happens with shifting to account for displacement
    dp_shift = 0
    dps = [-1,1]*dp_shift



plotopt = 6+1+1+1
plotopt = 7
;plotopt = 8
;plotopt = 9

;plotopt = 10
;plotopt = 11
plotopt = 7
;plotopt = 12
;plotopt = 13
plotopt = 14 ;2 col


zr1 = [-1,1]*3.5e20
zr1= [-1,1]*3e18*1.6


eps = 1*0

extag = ''
;extag = '_3'

if eps eq 1 then begin
    

	fname='~/temp.eps'
    xsize = 3.3
    ysize = 8;+2+2
    
    if plotopt eq 7 then begin
        fname = '~/m3dc1_alne_'+strtrim(shot,2)+'_2s_'+mtag+extag+'.eps'
        xsize = 3.3 & ysize =5+4
    endif		


    if plotopt eq 10 then begin
        fname = '~/m3dc1_thom_'+strtrim(shot,2)+'_avg_sym_line_'+mtag+extag+'.eps'
        xsize = 3.3 & ysize =4
    endif		

    if plotopt eq 11 then begin
        fname = '~/m3dc1_alne_refl_'+strtrim(shot,2)+'_'+mtag+extag+'.eps'
        xsize = 3.3 & ysize =5
    endif		

    if plotopt eq 12 then begin
        fname = '~/m3dc1_field_tor_'+strtrim(shot,2)+'_'+mtag+extag+'.eps'
        xsize = 8 & ysize =3.5
    endif		


    if plotopt eq 13 then begin
        fname = '~/m3dc1_refl_check_'+strtrim(shot,2)+'_'+mtag+extag+'.eps'
        xsize = 3.3 & ysize =3
    endif		


    if plotopt eq 14 then begin
        fname = '~/m3dc1_al_2c0ol_'+strtrim(shot,2)+'_'+mtag+extag+'.eps'
        xsize = 3.38 & ysize =4
    endif		

    		
	print,fname
    set_plot,'ps'
    device,/encapsulated,/color,/bold,xsize=xsize,set_font='Helvetica',$
    	/TT_FONT,ysize=ysize,font_size=6,preview=0,bits_per_pixel=14,$
		filename=fname,/inches;,xoffset=.5,yoffset=.5;/inches
	
	!p.font=0 & axisthick=!x.thick & !x.thick=2 & !y.thick=2
	!p.color=1 &!p.thick=2  & !x.ticklen=0.05

endif



if plotopt eq 1 then begin

    ;color_setup,1,/reverse
    col1 = blue
    col2 = red

    !P.multi=[0,1,3]

    plot,rv,cut_eq,xr = [2.2,2.3],linestyle=2,ytitle='Pe (Pa)'
    oplot,rv,cut_eq+cut_p1,col=col1
    oplot,rv,cut_eq+cut_p2,col=col2

    plot,rv,cut_p1,xr = !x.crange, yr = [-1,1]*250,/nodata,ytitle = 'Pe perturbation (Pa)'
    oplot,rv,cut_p1,col=col1
    oplot,rv,cut_p2,col=col2

    plot,rv,deriv(rv,cut_eq)/1e3,xr=!x.crange,linestyle=2,ytitle='Pe gradient (kPa/m)'
    oplot,rv,deriv(rv,cut_eq+cut_p1)/1e3,col=col1
    oplot,rv,deriv(rv,cut_eq+cut_p2)/1e3,col=col2

endif

    @greek_chars.pro
    @colors_kc.pro
    !p.background = white
    !p.color = black
    xr = [2.2,2.3]
    cols = [blue,red]
    cols = [red,blue]



if plotopt eq 2 then begin

    !p.multi=[0,2,2]
    for i = 0,nfields-1 do begin
        
        ;yr = [min(pertsl[*,i,*]),max(pertsl[*,i,*])]
        yr = [min(eqsl[*,i]),max(eqsl[*,i])]*2
    
        plot,rv,eqsl[*,i],title='a/'+labels[i],linestyle=2,yr=yr,xr=xr,xstyle=1
        oplot,[1,1]*rv_sep,[-1,1]*1e10,linestyle=1
        
        for j =0,nphis-1 do oplot,rv,pertsl[*,i,j],col=cols[j]
        
    endfor


endif

if plotopt eq 3 then begin

   !p.multi=[0,2,2]
    for i = 0,nfields-1 do begin
        
        ;yr = [min(pertsl[*,i,*]),max(pertsl[*,i,*])]
    
        plot,rv,eqs[*,i],title=fields[i],linestyle=2,xr=xr,xstyle=1
        
        for j =0,nphis-1 do oplot,rv,eqs[*,i]+perts[*,i,j],col=cols[j]
        
    endfor



endif

if plotopt eq 4 then begin

    !p.multi=[0,2,2]
    for i = 0,nfields-1 do begin
        
        yr = [min(perts[*,i,*]),max(perts[*,i,*])]
    
        plot,rv,eqs[*,i],title=fields[i],linestyle=2,xr=xr,xstyle=1,/nodata,yr=yr
        
        for j =0,nphis-1 do oplot,rv,perts[*,i,j],col=cols[j]
        
    endfor



endif


if plotopt eq 5 then begin

    !p.multi=[0,2,2]
    for i = 0,nfields-1 do begin
        
        yr = [min(perts[*,i,*]),max(perts[*,i,*])]
    
        plot,rv,deriv(rv,eqs[*,i]),title=fields[i],linestyle=2,xr=xr,xstyle=1
        
        for j =0,nphis-1 do oplot,rv,deriv(rv,eqs[*,i]+perts[*,i,j]),col=cols[j]
        
    endfor



endif

if plotopt eq 6 then begin

    rad = psinrv & xr = [0.8,1.05] & sepval = 1 & xtitle = psinstr
    ;rad = rv & xr = [2.1,2.3] & sepval = rv_sep & xtitle = 'R (cm)'
    i = where(fields eq 'den')
    ;i = where(fields eq 'ti')
    ;i = where(fields eq 'te')
    
    nplots=3
    xpl = 0.19 & xpml = 0.95
    extop=0.03 
    sep=0.01/2
    bas=0.1
    ytop=0.98-extop
    rang=ytop-bas
    text=rang/(nplots)    ;nis is # of plots
    pext=text-sep
    ypl=findgen(nplots)/(nplots-1)*(rang-text)+bas
    ypu=ypl+text-sep
    ;!y.ticklen = 0.005


    nplt=nplots-1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    plot,rad,eqs[*,i],linestyle=2,xr=xr,xstyle=5,ystyle=5,/nodata
    pow = floor(alog10(max(abs(!y.crange))))
    plot,rad,eqs[*,i]/(10.^pow),linestyle=2,xr=xr,xstyle=1,col=black,$
        xtickname=replicate(' ',6)
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.85),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N)'
    for j =0,nphis-1 do oplot,rad,(eqs[*,i]+perts[*,i,j])/(10.^pow),col=cols[j]

    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]
    plot,rad,deriv(rv,eqs[*,i]),linestyle=2,xr=xr,xstyle=5,ystyle=5,/nodata,/noerase
    pow = floor(alog10(max(abs(!y.crange))))
    plot,rad,deriv(rv,eqs[*,i])/(10.^pow),linestyle=2,xr=xr,xstyle=1,$
        /noerase,xtickname=replicate(' ',6)
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.85),gradstr+labels2[i]+' (x10!U'+strtrim(pow,2)+'!N)'
    for j =0,nphis-1 do oplot,rad,deriv(rv,eqs[*,i]+perts[*,i,j])/(10.^pow),col=cols[j]
    

    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]
    yr = [min(eqsl[*,i]),max(eqsl[*,i])]
    plot,rad,eqsl[*,i],linestyle=2,yr=yr,xr=xr,xstyle=1,$
        xtitle=xtitle,/noerase
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.85),'a/'+labels[i]
    for j =0,nphis-1 do oplot,rad,pertsl[*,i,j],col=cols[j]

  

endif

if plotopt eq 7 then begin

    rad = psinrv & xr = [0.85,1.05] & sepval = 1 & xtitle = psinstr
    ;rad = rv & xr = [2.1,2.3] & sepval = rv_sep & xtitle = 'R (cm)'
    i = where(fields eq 'den')
    ;i = where(fields eq 'ti')
    ;i = where(fields eq 'te')
    
    nplots=4+1+1
    ;xpl = 0.19 & xpml = 0.95-0.04
    xpl = 0.19-0.05 & xpml = 0.95-0.08
    extop=0.03 
    sep=0.01/2
    bas=0.1
    ytop=0.98-extop
    rang=ytop-bas
    text=rang/(nplots)    ;nis is # of plots
    pext=text-sep
    ypl=findgen(nplots)/(nplots-1)*(rang-text)+bas
    ypu=ypl+text-sep

    nplt=nplots

    ;------------------------------------
    ;   m3dc1 density perturbation in phi
    load_difct
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den') & pow = 3
    val = reform(perts_tot[*,i,*])
    vald = val*0
    for j = 0, ntot-1 do begin
        vald[*,j] = deriv(rv,val[*,j]);+eqs[*,i])
    endfor
    pixplot,val,psinrv,phi_machine_tot,xr=xr,zr=zr1,$
        ytitle=phistr+ ' (deg.)',xtickname=replicate(' ',6)
    nlev = 8
    levels = (findgen(nlev)/(nlev-1)*2-1)*zr1[1]
    ;contour,vald,psinrv,phi_machine_tot,xr=xr,zr=[-1,1]*3.5e20,$
    ;    ytitle=phistr+ ' (deg.)',xtickname=replicate(' ',6),/fill,levels=levels
    ;contour,vald,psinrv,phi_machine_tot,/over,levels=levels
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    oplot,besstr.psf_avgs,besstr.phipos*!radeg,col=black
    oplot,besstr.psf_avgs,besstr.phipos*!radeg,col=black,psym=sym(1)

    psfi = 3
    bpsin2d = quickpsinorm(besstr.g)
    int1 = interpol(lindgen(besstr.g.mw),besstr.g.r,besstr.psfs[psfi].r/1e2)
    int2 = interpol(lindgen(besstr.g.mh),besstr.g.z,0)
    psfpsin = interpolate(bpsin2d,int1,int2,/grid)
    oplot,psfpsin,besstr.psfs[psfi].psf_r*360
    ;oplot,psinrv,interpol(besstr.psfs[4].psf_r,besstr.psfs[4].r/1e2,rv)*360

	;side color bar
	pow = 18
	zr1 = zr1/10.^pow
	!p.position[0]=0.97 & !p.position[2]=0.99
	ncbp = 256
    pixplot,transpose(rebin(zr1,ncbp,2)),[0,1],rebin(zr1,ncbp),/noerase,$
    	xticklen=0.0001,xtickname=replicate(' ',6),yticklen=.2,$
		ytickinterval=2,charsize=1,$
		ytitle=deltastr+labels2[i]+ ' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)
		;ytitle=gradstr+deltastr+labels2[i]+ ' (x10!U'+strtrim(pow,2)+$
		;'!N m!U-4!N)',



    @colors_kc.pro
    !p.background = white
    !p.color = black
    ;------------------------------------
    ;   m3dc1 te profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'te') & pow = 3
    plot,rad,eqs[*,i]/(10.^pow),linestyle=2,xr=xr,xstyle=1,col=black,/noerase,$
        xtickname=replicate(' ',6),yr = [0,1.5],ytickinterval = 1,ystyle=1
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    ;xyouts,xcr(0.05),ycr(0.85),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N)'
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (keV)'    
    for j =0,nphis-1 do oplot,dps[j]+rad,(eqs[*,i]+perts[*,i,j])/(10.^pow),col=cols[j]
    xyouts,xcr(0.05),ycr(0.05),tag,color=black,charsize=0.5

    ;oplot,p.psi.tanh.psi_te,p.psi.tanh.t_e,col=green

    ;------------------------------------
    ;   m3dc1 ne profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den') & pow = 19
    plot,rad,eqs[*,i]/(10.^pow),linestyle=2,xr=xr,xstyle=1,col=black,$
        xtickname=replicate(' ',6),yr = [0,3],ytickinterval = 1,/noerase
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'
    for j =0,nphis-1 do oplot,dps[j]+rad,(eqs[*,i]+perts[*,i,j])/(10.^pow),col=cols[j]

    ;------------------------------------
    ;   m3dc1 a/Lne profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den')
    plot,rad,ga_smooth(abs(eqsl[*,i]),rad,0.005),linestyle=2,xr=xr,xstyle=1,$
        /noerase,ytickinterval=50,xtickname=replicate(' ',6)
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.85),'a/'+labels[i]
    for j =0,nphis-1 do oplot,dps[j]+rad,ga_smooth(abs(pertsl[*,i,j]),rad,0.005,/boxcar),col=cols[j]
    ;for j =0,nphis-1 do oplot,rad,abs(pertsl[*,i,j]),col=cols[j]

    ;------------------------------------
    ;   m3dc1 a/Lti profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'ti')
    plot,rad,ga_smooth(abs(eqsl[*,i]),rad,0.001),linestyle=2,xr=xr,xstyle=1,$
        xtickname=replicate(' ',6),/noerase,ytickinterval=50
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.85),'a/'+labels[i]
    ;for j =0,nphis-1 do oplot,rad,abs(pertsl[*,i,j]),col=cols[j]
    for j =0,nphis-1 do oplot,dps[j]+rad,ga_smooth(abs(pertsl[*,i,j]),rad,0.001,/boxcar),col=cols[j]

    ;------------------------------------
    ;  BES dI/I for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i0 = where(besstr.tbase eq time)
    i1 = i0+1
    besis = [i0,i1]
    syms = [5,1]
    yerr= replicate(0,6)
    plot,[0,0],/nodata,xr=xr,yr=[0,0.7],ystyle=1,/noerase,xtitle=psinstr,ytickinterval=0.5
    xyouts,xcr(0.05),ycr(0.85),deltastr+'I/I (%)'
    for j=0,n_elements(besis)-1 do begin
        ii = besis[j]
        oploterror,besstr.psf_avgs,besstr.an_tot[*,ii],besstr.psf_xerrs,yerr,$
                    errcol=cols[j],col=cols[j]
        oplot,besstr.psf_avgs,besstr.an_tot[*,ii],psym=sym(syms[j]),col=cols[j]
        oplot,besstr.psf_avgs,besstr.an_tot[*,ii],psym=sym(syms[j]+5),col=black
        xyouts,xcr(0.05+0.2*j),ycr(0.7),strtrim(round(besstr.tbase[ii]),2),$
            col=cols[j],charsize=1
    endfor
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    oplot,psfpsin,besstr.psfs[psfi].psf_r*!y.crange[1]



endif


if plotopt eq 8 then begin

    !p.multi=[0,4,1]
    ;load_difct
    load_mod33ct
    xr=[0.8,1.05]
    
    !p.charsize=2.5
    
    for i = 0, nfields -1 do begin
        
        val = reform(perts_tot[*,i,*])
        vald = val*0
        for j = 0, ntot-1 do begin
            vald[*,j] = deriv(rv,val[*,j]);+eqs[*,i])
        endfor
        ;pixplotcb,reform(perts_tot[*,i,*]),psinrv,phi_machine_tot,xr=xr,$
        ;    cbarlabel = deltastr+labels2[i],/autoz
        pixplotcb,vald,psinrv,phi_machine_tot,xr=xr,$
            cbarlabel = gradstr+deltastr+labels2[i],/autoz,xtitle=psinstr,$
            ytitle=phistr+ ' (deg.)'
        ;pixplotcb,abs(reform(pertsl_tot[*,i,*])),psinrv,phi_machine_tot,xr=xr,$
        ;    cbarlabel = 'a/'+labels[i],/autoz
        oplot,[0,10],[1,1]*146
        oplot,[1,1],[0,400],linestyle=1
    
    endfor
    

endif

if plotopt eq 9 then begin

    @colors_kc.pro
    !p.background = white
    !p.color = black

    dsym = 3
!p.multi=[0,2,2]
    
    i = where(fields eq 'den') & pow = 19
    plot,zt,eqst[*,i]/(10.^pow),linestyle=2,xstyle=1,col=black,$
        yr = [0,3],ytickinterval = 1,xr=[0.6,0.77],xtitle='Z (m)';,xtitle=psin
    oplot,[1,1]*zt_sep,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'
    for j =0,nphist-1 do oplot,zt,(eqst[*,i]+pertst[*,i,j])/(10.^pow),col=cols[j]

;    oplot,zt,p1ne,col=purple,linestyle=2,thick=3
;    oplot,zt,p2ne,col=ltblue,linestyle=2,thick=3

    oplot,p1th_z,p1.psi.dat.n_e*10,col=red,psym=dsym
    oplot,p2th_z,p2.psi.dat.n_e*10,col=blue,psym=dsym
   

    for j=0,ntags-1 do begin
        
        r=EXECUTE('pj = profiles.'+tags[j])
        ;oplot,pj.z/1e2,pj.dens,psym=3,col=cols[j]
        ;oplot,pj.z_dens_fit/1e2,pj.dens_fit_on_z,col=cols[j],thick=3,linestyle=3
        
    endfor

    plot,psinzt,eqst[*,i]/(10.^pow),linestyle=2,xstyle=1,col=black,$
        yr = [0,3],ytickinterval = 1,xr=[0.85,1.05],xtitle=psinstr;,xtitle=psin
    oplot,[1,1],[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'
    for j =0,nphist-1 do oplot,psinzt,(eqst[*,i]+pertst[*,i,j])/(10.^pow),col=cols[j]

    ;oplot,p1psin,p1ne,col=purple,thick=3,linestyle=2
    ;oplot,p2psin,p2ne,col=ltblue,thick=3,linestyle=2
    oplot,p1.psi.tanh.psi_ne-qres[gqr[0]].adjust_zts+exshift,p1.psi.tanh.n_e*10,col=purple,linestyle=2,thick=3
    oplot,p2.psi.tanh.psi_ne-qres[gqr[1]].adjust_zts+exshift,p2.psi.tanh.n_e*10,col=ltblue,linestyle=2,thick=3
;    oplot,p1.psi.tanh.psi_ne,p1.psi.tanh.n_e*10,col=purple,linestyle=2,thick=3
;    oplot,p2.psi.tanh.psi_ne,p2.psi.tanh.n_e*10,col=ltblue,linestyle=2,thick=3

    i = where(fields eq 'te') & pow = 3
    plot,zt,eqst[*,i]/(10.^pow),linestyle=2,xstyle=1,col=black,$
        yr = [0,2],ytickinterval = 1,xr=[0.6,0.77],xtitle='Z (m)';,xtitle=psin
    oplot,[1,1]*zt_sep,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'
    for j =0,nphist-1 do oplot,zt,(eqst[*,i]+pertst[*,i,j])/(10.^pow),col=cols[j]

    for j=0,ntags-1 do begin
        
        r=EXECUTE('pj = profiles.'+tags[j])
        ;oplot,pj.z/1e2,pj.temp,psym=3,col=cols[j]
        ;oplot,pj.z_temp_fit/1e2,pj.temp_fit_on_z,col=cols[j],thick=3,linestyle=3
    
    endfor

;    oplot,zt,p1te,col=purple,linestyle=2,thick=3
;    oplot,zt,p2te,col=ltblue,linestyle=2,thick=3

    oplot,p1th_z,p1.psi.dat.t_e,col=red,psym=dsym
    oplot,p2th_z,p2.psi.dat.t_e,col=blue,psym=dsym
    

    plot,psinzt,eqst[*,i]/(10.^pow),linestyle=2,xstyle=1,col=black,$
        yr = [0,2],ytickinterval = 1,xr=[0.85,1.05],xtitle=psinstr;,xtitle=psin
    oplot,[1,1],[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'
    for j =0,nphist-1 do oplot,psinzt,(eqst[*,i]+pertst[*,i,j])/(10.^pow),col=cols[j]
    
    ;oplot,p1psin,p1te,col=purple,thick=3,linestyle=2
    ;oplot,p1psin,p2te,col=ltblue,thick=3,linestyle=2
    oplot,p1.psi.tanh.psi_te-qres[gqr[0]].adjust_zts+exshift,p1.psi.tanh.t_e,col=purple,linestyle=2,thick=3
    oplot,p2.psi.tanh.psi_te-qres[gqr[1]].adjust_zts+exshift,p2.psi.tanh.t_e,col=ltblue,linestyle=2,thick=3
;    oplot,p1.psi.tanh.psi_te,p1.psi.tanh.t_e,col=purple,linestyle=2,thick=3
;    oplot,p2.psi.tanh.psi_te,p2.psi.tanh.t_e,col=ltblue,linestyle=2,thick=3
    
endif

if plotopt eq 10 then begin


    cms = 1e2 & xtici = 1 & xtit = 'Z (cm)'
    ;cms = 1 & xtici = 0.05 & xtit = 'Z (m)'
    
    xr = [0.65,0.77]
    xr = [0.68,0.77]
    xr = [0.69,0.76]
    xr = [0.70,0.75]*cms
    

    dsym=sym(1+5)
    dsym2 = sym(5+5)
    @colors_kc.pro
    !p.background = white
    !p.color = black

    nplots=2
    xpl = 0.19-0.05 & xpml = 0.95
    extop=0.03 
    sep=0.01/2
    bas=0.1
    ytop=0.98-extop
    rang=ytop-bas
    text=rang/(nplots)    ;nis is # of plots
    pext=text-sep
    ypl=findgen(nplots)/(nplots-1)*(rang-text)+bas
    ypu=ypl+text-sep

    nplt=nplots

    !x.ticklen=0.02
    ;------------------------------------
    ;   m3dc1 density v data
    yr = [0,2.5]
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den') & pow = 19
    plot,zt*cms,eqst[*,i]/(10.^pow),linestyle=2,xstyle=1,col=black,$
        yr = yr,ytickinterval = 1,xr=xr,xtickname=replicate(' ',6),$
        xtickinterval=xtici,ystyle=1
    oplot,[1,1]*zt_sep*cms,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'
    for j =0,nphist-1 do oplot,zt*cms,(eqst[*,i]+pertst[*,i,j])/(10.^pow),col=cols[j]

    oplot,zt*cms,p1ne,col=purple,linestyle=2,thick=3
    oplot,zt*cms,p2ne,col=ltblue,linestyle=2,thick=3

;    oplot,p1th_z*cms,p1.psi.dat.n_e*10,col=red,psym=3
;    oplot,p2th_z*cms,p2.psi.dat.n_e*10,col=blue,psym=3

    oploterror,p1thzh*cms,p1neh,p1thzhe*cms,p1nehe,errcol=red,col=red,psym=3,/nohat
;    oplot,p1thzh*cms,p1neh,psym=dsym,col=red,symsize=0.75
    oploterror,p2thzh*cms,p2neh,p2thzhe*cms,p2nehe,errcol=blue,col=blue,psym=3,/nohat
;    oplot,p2thzh*cms,p2neh,psym=dsym2,col=blue,symsize=0.75

    
    ;------------------------------------
    ;   m3dc1 te v data
    yr =[0,2]
    yr =[0,1.2]
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'te') & pow = 3
    plot,zt*cms,eqst[*,i]/(10.^pow),linestyle=2,xstyle=1,col=black,/noerase,$
        yr = yr,ytickinterval = 1,xr=xr,xtitle=xtit,ystyle=1,$
        xtickinterval=xtici;,xtitle=psin
    oplot,[1,1]*zt_sep*cms,[-1,1]*1e10,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.4),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'
    for j =0,nphist-1 do oplot,zt*cms,(eqst[*,i]+pertst[*,i,j])/(10.^pow),col=cols[j]

    oplot,zt*cms,p1te,col=purple,linestyle=2,thick=3
    oplot,zt*cms,p2te,col=ltblue,linestyle=2,thick=3

;    oplot,p1th_z*cms,p1.psi.dat.t_e,col=red,psym=3
;    oplot,p2th_z*cms,p2.psi.dat.t_e,col=blue,psym=3

;    oplot,p1thzh*cms,p1teh,psym=dsym,col=red
;    errplot,p1thzh*cms,p1teh-p1tehe,p1teh+p1tehe,col=red
;    oplot,p2thzh*cms,p2teh,psym=dsym,col=blue
;    errplot,p2thzh*cms,p2teh-p2tehe,p2teh+p2tehe,col=blue
 
    oploterror,p1thzh*cms,p1teh,p1thzhe*cms,p1tehe,errcol=red,col=red,psym=3,/nohat
;    oplot,p1thzh*cms,p1teh,psym=dsym,col=red,symsize=0.75
    oploterror,p2thzh*cms,p2teh,p2thzhe*cms,p2tehe,errcol=blue,col=blue,psym=3,/nohat
;    oplot,p2thzh*cms,p2teh,psym=dsym2,col=blue,symsize=0.75
  

endif

if plotopt eq 11 then begin

    rad = psinrv & xr = [0.85,1.05] & sepval = 1 & xtitle = psinstr

    nplots=2+1+1
    xpl = 0.19-0.05 & xpml = 0.95
    xpl = 0.19-0.05 & xpml = 0.95-0.08

    extop=0.03-0.01 
    sep=0.01/2
    bas=0.1
    ytop=0.98-extop
    rang=ytop-bas
    text=rang/(nplots)    ;nis is # of plots
    pext=text-sep
    ypl=findgen(nplots)/(nplots-1)*(rang-text)+bas
    ypu=ypl+text-sep

    nplt=nplots

    ;------------------------------------
    ;   m3dc1 density perturbation in phi
    load_difct
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den') & pow = 3
    val = reform(perts_tot[*,i,*])
    vald = val*0
    for j = 0, ntot-1 do begin
        vald[*,j] = deriv(rv,val[*,j]);+eqs[*,i])
    endfor
    ;zr1 = [-1,1]*3.5e20
    ;zr1= [-1,1]*3e18
    pixplot,val,psinrv,phi_machine_tot,xr=xr,zr=zr1,$
        ytitle=phistr+ ' (deg.)',xtickname=replicate(' ',6)
    nlev = 8
    levels = (findgen(nlev)/(nlev-1)*2-1)*zr1[1]
    ;contour,vald,psinrv,phi_machine_tot,xr=xr,zr=[-1,1]*3.5e20,$
    ;    ytitle=phistr+ ' (deg.)',xtickname=replicate(' ',6),/fill,levels=levels
    ;contour,vald,psinrv,phi_machine_tot,/over,levels=levels
    oplot,[0,10],[1,1]*255,linestyle=2
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black

	;side color bar
	pow = 18
	zr1 = zr1/10.^pow
	!p.position[0]=0.97 & !p.position[2]=0.99
	ncbp = 256
    pixplot,transpose(rebin(zr1,ncbp,2)),[0,1],rebin(zr1,ncbp),/noerase,$
    	xticklen=0.0001,xtickname=replicate(' ',6),yticklen=.2,$
		ytickinterval=2,charsize=1,$
		ytitle=deltastr+labels2[i]+ ' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)
		;ytitle=gradstr+deltastr+labels2[i]+ ' (x10!U'+strtrim(pow,2)+$
		;'!N m!U-4!N)',

    @colors_kc.pro
    !p.background = white
    !p.color = black
    ;------------------------------------
    ;   m3dc1 te profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'te') & pow = 3
    plot,rad,eqs[*,i]/(10.^pow),linestyle=2,xr=xr,xstyle=1,col=black,/noerase,$
        xtickname=replicate(' ',6),yr = [0,1.5],ytickinterval = 1,ystyle=1
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    ;xyouts,xcr(0.05),ycr(0.85),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N)'
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (keV)'    
    for j =0,nphis-1 do oplot,dps[j]+rad,(eqs[*,i]+pertsr[*,i,j])/(10.^pow),col=cols[j]
    xyouts,xcr(0.05),ycr(0.05),tag,color=black,charsize=0.5

    ;oplot,p.psi.tanh.psi_te,p.psi.tanh.t_e,col=green

    ;------------------------------------
    ;   m3dc1 ne profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den') & pow = 19
    plot,rad,eqs[*,i]/(10.^pow),linestyle=2,xr=xr,xstyle=1,col=black,$
        xtickname=replicate(' ',6),yr = [0,3],ytickinterval = 1,/noerase
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'
    for j =0,nphis-1 do oplot,dps[j]+rad,(eqs[*,i]+pertsr[*,i,j])/(10.^pow),col=cols[j]


    ;------------------------------------
    ;   m3dc1 a/Lne profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den')
    plot,rad,ga_smooth(abs(eqsl[*,i]),rad,0.005),linestyle=2,xr=xr,xstyle=1,$
        /noerase,ytickinterval=50,xtitle=psinstr
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.85),'a/'+labels[i]
    for j =0,nphis-1 do oplot,dps[j]+rad,ga_smooth(abs(pertslr[*,i,j]),rad,0.005,/boxcar),col=cols[j]
    ;for j =0,nphis-1 do oplot,rad,abs(pertsl[*,i,j]),col=cols[j]
    for j =0,nphis-1 do xyouts,xcr(0.05),ycr(0.6-j*0.2),phistr+' = '+$
        strtrim(360+(30-phisr[j]),2),col=cols[j]


endif

if plotopt eq 12 then begin

    rad = psinrv & xr = [0.85,1.05] & sepval = 1 & xtitle = psinstr

    !p.multi=[0,4,1]
    !p.charsize = 3
    if !d.name eq 'PS' then !p.charsize=1.5
   ;------------------------------------
    ;   m3dc1 density perturbation in phi
    load_difct
    ;nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den') & pow = 3
    
    for i = 0,nfields-1 do begin
    
    val = reform(perts_tot[*,i,*])
    vald = val*0
    for j = 0, ntot-1 do begin
        vald[*,j] = deriv(rv,val[*,j]);+eqs[*,i])
    endfor
    zr1 = [-1,1]*3.5e20
    zr1= [-1,1]*3e18
    pixplotcb,val,psinrv,phi_machine_tot,xr=xr,/autoz,$
        ytitle=phistr+ ' (deg.)',xtitle=xtitle,cbarlabel=fields[i]
    ;nlev = 8
    ;levels = (findgen(nlev)/(nlev-1)*2-1)*zr1[1]
    ;contour,vald,psinrv,phi_machine_tot,xr=xr,zr=[-1,1]*3.5e20,$
    ;    ytitle=phistr+ ' (deg.)',xtickname=replicate(' ',6),/fill,levels=levels
    ;contour,vald,psinrv,phi_machine_tot,/over,levels=levels
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    oplot,besstr.psf_avgs,besstr.phipos*!radeg,col=black
    oplot,besstr.psf_avgs,besstr.phipos*!radeg,col=black,psym=sym(1)

    oplot,[0,10],[1,1]*255,linestyle=2,thick=3

    endfor

    xyouts,0.2,0.01,file,/normal

    !p.multi=[0,4,1]
    !p.charsize = 1
endif


if plotopt eq 13 then begin

    colst = [black,red]
    colst = cols
    i = where(fields eq 'den') & pow = 19

    ;plot,[0,0],/nodata,xr=[2.25,2.3],yr=[0,2.5e19]
    ;oploterror,pr1r,pr1d,pr1re,pr1de,errcol=red,col=red,psym=3,/nohat
    ;oploterror,pr2r,pr2d,pr2re,pr2de,errcol=blue,col=blue,psym=3,/nohat
    ;for j =0,nphisr-1 do oplot,rv,(eqs[*,i]+pertsr[*,i,j])/(10.^pow),col=cols[j]

    plot,[0,0],/nodata,xr=[0.85,1.05],xstyle=1,yr=[0,3],xtitle=psinstr
    oploterror,pr1psin,pr1d/1e19,pr1psine,pr1de/1e19,errcol=colst[0],col=colst[0],psym=3,/nohat
    oploterror,pr2psin,pr2d/1e19,pr2psine,pr2de/1e19,errcol=colst[1],col=colst[1],psym=3,/nohat
    for j =0,nphisr-1 do oplot,psinrv,(eqs[*,i]+pertsr[*,i,j])/(10.^pow),col=colst[j],linestyle=2
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'

;oplot,rv,pertsr[*,i,0]+eqs[*,1] 

endif

if plotopt eq 14 then begin
; two column radial blah blah blah

    rad = psinrv & xr = [0.85,1.05] & sepval = 1 & xtitle = psinstr
    ;rad = rv & xr = [2.1,2.3] & sepval = rv_sep & xtitle = 'R (cm)'

    ;columns!
    nplots=4 ;per column

    x_sep = 0.09
    xpl_1 = 0.05 & xpml_2 = 0.915
    wid = (xpml_2 - xpl_1 - x_sep)/2
    xpml_1 = (xpl_1 + wid) & xpl_2 = (xpml_1 + x_sep)
    extop=0.00
    sep=0.01/2
    bas=0.1
    ytop=0.98-extop
    rang=ytop-bas
    text=rang/(nplots)    ;nis is # of plots
    pext=text-sep
    ypl=findgen(nplots)/(nplots-1)*(rang-text)+bas
    ypu=ypl+text-sep

  
    ;--- col 1 ------------------------------------------------
    xpl = xpl_1 & xpml = xpml_1
    nplt=nplots
    ;----------------------------------------------------------
    
    @colors_kc.pro
    !p.background = white
    !p.color = black
    ;------------------------------------
    ;   m3dc1 p profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'p') & pow = 3
    plot,rad,eqs[*,i]/(10.^pow),linestyle=2,xr=xr,xstyle=1,col=black,$
        xtickname=replicate(' ',6),ytickinterval = 5,ystyle=1
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    ;xyouts,xcr(0.05),ycr(0.85),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N)'
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (kPa)'    
    for j =0,nphis-1 do oplot,dps[j]+rad,(eqs[*,i]+perts[*,i,j])/(10.^pow),col=cols[j]
    xyouts,xcr(0.05),ycr(0.05),tag,color=black,charsize=0.5

    ;------------------------------------
    ;   m3dc1 te profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'te') & pow = 3
    plot,rad,eqs[*,i]/(10.^pow),linestyle=2,xr=xr,xstyle=1,col=black,$
        xtickname=replicate(' ',6),yr = [0,1.5],ytickinterval = 1,ystyle=1,/noerase
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    ;xyouts,xcr(0.05),ycr(0.85),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N)'
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (keV)'    
    for j =0,nphis-1 do oplot,dps[j]+rad,(eqs[*,i]+perts[*,i,j])/(10.^pow),col=cols[j]
    xyouts,xcr(0.05),ycr(0.05),tag,color=black,charsize=0.5


;    ;------------------------------------
;    ;   m3dc1 a/LTe profile for two phases
;    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
;    i = where(fields eq 'te')
;    plot,rad,ga_smooth(abs(eqsl[*,i]),rad,0.005),linestyle=2,xr=xr,xstyle=1,$
;        /noerase,ytickinterval=100,xtickname=replicate(' ',6)
;    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
;    xyouts,xcr(0.05),ycr(0.85),'a/'+labels[i]
;    for j =0,nphis-1 do oplot,dps[j]+rad,ga_smooth(abs(pertsl[*,i,j]),rad,0.005,/boxcar),col=cols[j]
;    ;for j =0,nphis-1 do oplot,rad,abs(pertsl[*,i,j]),col=cols[j] 

    ;------------------------------------
    ;   m3dc1 ti profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'ti') & pow = 3
    plot,rad,eqs[*,i]/(10.^pow),linestyle=2,xr=xr,xstyle=1,col=black,$
        xtickname=replicate(' ',6),yr = [0,2],ytickinterval = 1,ystyle=1,/noerase
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    ;xyouts,xcr(0.05),ycr(0.85),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N)'
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (keV)'    
    for j =0,nphis-1 do oplot,dps[j]+rad,(eqs[*,i]+perts[*,i,j])/(10.^pow),col=cols[j]
    xyouts,xcr(0.05),ycr(0.05),tag,color=black,charsize=0.5

    ;------------------------------------
    ;   m3dc1 a/LTi profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'ti')
    plot,rad,ga_smooth(abs(eqsl[*,i]),rad,0.005),linestyle=2,xr=xr,xstyle=1,$
        /noerase,ytickinterval=50,xtitle=psinstr
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.85),'a/'+labels[i]
    for j =0,nphis-1 do oplot,dps[j]+rad,ga_smooth(abs(pertsl[*,i,j]),rad,0.005,/boxcar),col=cols[j]
    ;for j =0,nphis-1 do oplot,rad,abs(pertsl[*,i,j]),col=cols[j]


    ;--- col 2 ------------------------------------------------
    xpl = xpl_2 & xpml = xpml_2
    nplt=nplots
    ;----------------------------------------------------------

    ;------------------------------------
    ;   m3dc1 density perturbation in phi
    load_difct
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den') & pow = 3
    val = reform(perts_tot[*,i,*])
    vald = val*0
    for j = 0, ntot-1 do begin
        vald[*,j] = deriv(rv,val[*,j]);+eqs[*,i])
    endfor
    pixplot,val,psinrv,phi_machine_tot,xr=xr,zr=zr1,ytickinterval=90,yminor=9,$
        ytitle=phistr+ ' (deg.)',xtickname=replicate(' ',6),/noerase
    nlev = 8
    levels = (findgen(nlev)/(nlev-1)*2-1)*zr1[1]
    ;contour,vald,psinrv,phi_machine_tot,xr=xr,zr=[-1,1]*3.5e20,$
    ;    ytitle=phistr+ ' (deg.)',xtickname=replicate(' ',6),/fill,levels=levels
    ;contour,vald,psinrv,phi_machine_tot,/over,levels=levels
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    oplot,besstr.psf_avgs,besstr.phipos*!radeg,col=black
    oplot,besstr.psf_avgs,besstr.phipos*!radeg,col=black,psym=sym(1)

    psfi = 3
    bpsin2d = quickpsinorm(besstr.g)
    int1 = interpol(lindgen(besstr.g.mw),besstr.g.r,besstr.psfs[psfi].r/1e2)
    int2 = interpol(lindgen(besstr.g.mh),besstr.g.z,0)
    psfpsin = interpolate(bpsin2d,int1,int2,/grid)
    oplot,psfpsin,besstr.psfs[psfi].psf_r*360
    ;oplot,psinrv,interpol(besstr.psfs[4].psf_r,besstr.psfs[4].r/1e2,rv)*360

	;side color bar
	pow = 18
	zr1 = zr1/10.^pow
	!p.position[0]=0.99 & !p.position[2]=1.
	ncbp = 256
    pixplot,transpose(rebin(zr1,ncbp,2)),[0,1],rebin(zr1,ncbp),/noerase,$
    	xticklen=0.0001,xtickname=replicate(' ',6),yticklen=.2,$
		ytickinterval=2,charsize=1,$
		ytitle=deltastr+labels2[i]+ ' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)
		;ytitle=gradstr+deltastr+labels2[i]+ ' (x10!U'+strtrim(pow,2)+$
		;'!N m!U-4!N)',


    @colors_kc.pro
    !p.background = white
    !p.color = black

    ;------------------------------------
    ;   m3dc1 ne profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den') & pow = 19
    plot,rad,eqs[*,i]/(10.^pow),linestyle=2,xr=xr,xstyle=1,col=black,$
        xtickname=replicate(' ',6),yr = [0,3],ytickinterval = 1,/noerase
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.6),labels2[i]+' (x10!U'+strtrim(pow,2)+'!N m!U-3!N)'
    for j =0,nphis-1 do oplot,dps[j]+rad,(eqs[*,i]+perts[*,i,j])/(10.^pow),col=cols[j]

    ;------------------------------------
    ;   m3dc1 a/Lne profile for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i = where(fields eq 'den')
    plot,rad,ga_smooth(abs(eqsl[*,i]),rad,0.005),linestyle=2,xr=xr,xstyle=1,$
        /noerase,ytickinterval=50,xtickname=replicate(' ',6)
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    xyouts,xcr(0.05),ycr(0.85),'a/'+labels[i]
    for j =0,nphis-1 do oplot,dps[j]+rad,ga_smooth(abs(pertsl[*,i,j]),rad,0.005,/boxcar),col=cols[j]
    ;for j =0,nphis-1 do oplot,rad,abs(pertsl[*,i,j]),col=cols[j]

    ;------------------------------------
    ;  BES dI/I for two phases
    nplt-=1 & !P.position=[xpl,ypl[nplt],xpml,ypu[nplt]]	
    i0 = where(besstr.tbase eq time)
    i1 = i0+1
    besis = [i0,i1]
    syms = [5,1]
    yerr= replicate(0,6)
    plot,[0,0],/nodata,xr=xr,yr=[0,0.7],ystyle=1,/noerase,xtitle=psinstr,ytickinterval=0.5
    xyouts,xcr(0.05),ycr(0.85),deltastr+'I/I (%)'
    for j=0,n_elements(besis)-1 do begin
        ii = besis[j]
        oploterror,besstr.psf_avgs,besstr.an_tot[*,ii],besstr.psf_xerrs,yerr,$
                    errcol=cols[j],col=cols[j]
        oplot,besstr.psf_avgs,besstr.an_tot[*,ii],psym=sym(syms[j]),col=cols[j]
        oplot,besstr.psf_avgs,besstr.an_tot[*,ii],psym=sym(syms[j]+5),col=black
        xyouts,xcr(0.05+0.2*j),ycr(0.7),strtrim(round(besstr.tbase[ii]),2),$
            col=cols[j],charsize=1
    endfor
    oplot,[1,1]*sepval,[-1,1]*1e30,linestyle=1,col=black
    oplot,psfpsin,besstr.psfs[psfi].psf_r*!y.crange[1]



endif



if eps eq 1 then begin

   device,/close
   set_plot,'x'
   !p.font=-1
   !x.thick=0
   !y.thick=0
   !p.thick=0
   !x.ticklen=0
   spawn,'evince '+fname+' &'
endif


!p.charsize=1
!p.multi=[0,1,1]
!p.position=[0,0,0,0]

;endfor
;endfor


end
