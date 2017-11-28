; Load an M3D-C1 output file and dump the relevant data into an IDL .sav file

save_fileloc = '/u/wilcoxr/idl/dens_xsection.sav'

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

linfac = -2.5*4/!pi
ntor=3
pts = 501

; Radial profile (Z=0, R range)
nprof = 701
zv = fltarr(nprof)
rv = findgen(nprof)/(nprof-1)*0.3+2.1

;rrange = [2,2.5]
;zrange = [-0.5,0.5]
;rrange = [1.9,2.5]
;zrange = [-0.2,0.9]
rrange = [1.05,2.5]
zrange = [-1.3,1.05]
nr = 309  ; 5 mm spacing
nz = 471

;r_locs = findgen(nr)

fields = ['den']
;fields = ['te','den','ti','pe','p']


nfields = n_elements(fields)

phi_peak = 30.
phi_in = (((30 - phi_peak) + 360) mod 360)




nzt = 301
ztr = [-6.57,82.58]/1e2
zt = findgen(nzt)/(nzt-1)*(ztr[1]-ztr[0])+ztr[0]
rt = replicate(1.94,nzt)

if getdat eq 1 then begin

    p = get_profdb(shot,time,prunid)
    p1 = get_profdb(shot,time,prunid)
    p2 = get_profdb(shot,3850,'bw01')
    
    eqs = fltarr(nprof,nfields)
    perts = fltarr(nprof,nfields)
    eqst = fltarr(nzt,nfields)


    for i = 0,nfields-1 do begin

        ;rese = read_field(fields[i],r,z,t,slice=eqslice,filename=file,points=pts,$
        ;                  /mks,/linear,/equilibrium,rrange=rrange,zrange=zrange)

        ;int1 = interpol(lindgen(pts),r,rv)
        ;int2 = interpol(lindgen(pts),z,zv)
        ;eqs[*,i] = interpolate(reform(rese),int1,int2)

        res = read_field(fields[i],r,z,t,slice=pslice,filename=file,points=pts,/mks,$
                         /linear,/complex,rrange=rrange,zrange=zrange)

        imag = imaginary(reform(res))
        real = real_part(reform(res))   ;same as resreal...
        amp = linfac*sqrt(imag^2+real^2)
        phase = atan(reform(res),/phase)

        
        int1 = interpol(lindgen(pts),r,rv)
        int2 = interpol(lindgen(pts),z,zv)
        
        
        perts[*,i] = interpolate(amp*cos(phase+phi_in/!radeg*ntor),int1,int2)
        
    endfor

endif


    
save,r,z,eqs,perts,eqst,fields,res,filename=save_fileloc


end
