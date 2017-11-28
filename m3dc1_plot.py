import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import scipy.io as si
import scipy.interpolate as scinter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import h5py


# plt.rcParams.update({'font.family': 'Arial'})
plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'figure.facecolor': 'w'})
plt.rcParams.update({'mathtext.default': 'regular'})

a_minor = 0.591  # minor radius in m, for calculating a/Lne
# R_LCFS = 2.2735

OMFIT_M3DC1_folder = '/fusion/projects/results/m3dc1_results/wilcoxr/data/OMFIT/DIII-D/'


#-----------------------------------------------------------------------

class M3DC1_xsection(object):
    def __init__(self, sav_fileloc = '/home/wilcoxr/m3dc1/169938/xsection_169938_4950.sav',
                 verbose = False):
    # def __init__(self, sav_fileloc = '/home/wilcoxr/m3dc1/profile_savs/170007/xsection_170007_3625.sav',
    #              verbose = False):
        """
        Requires previously saved IDL .sav file containing profile variables
        I generate this data with this (on venus): /u/wilcoxr/idl/load_m3dc1_single_n_fullxsection.pro

        Gets the imaginary fields so you can evaluate at any toroidal location and U-L phasing
        ** Only does it for 1 'fields' value for now

        """
        
        idlsave_dict = si.readsav(sav_fileloc, python_dict = True, verbose = verbose)
        
        # fields = ['te','den','ti','pe','p'] (usually ordered like this in arrays)
        
        
        self.fields = idlsave_dict['fields']  # fields that are in the IDL sav file ('te', 'ne', 'ti', 'p', 'pe'),
                                              # and which index corresponds to which
        self.r = idlsave_dict['r']  # grid locations (401)
        self.z = idlsave_dict['z']
        self.psin2d = np.transpose(np.reshape(idlsave_dict['psinrz'], (len(self.r), len(self.z))))  # psinorm on R, Z grid
        
        self.eqs = {}  # arrays are transposed from IDL->Python, need to flip them here
        self.perts = {}
        
        for f in self.fields:
            if f == 'den':
                f = 'ne'
            
            self.eqs[f] = np.transpose(idlsave_dict[f + '_eq'][:, :, 0])  # raw equilibrium values (should be real)
            
            for row in ['_u', '_l']:
                # self.perts[f + row] = idlsave_dict[f + '_pert_comp' + row][:, :, 0]
                self.perts[f + row] = np.transpose(idlsave_dict[f + '_pert_comp' + row][:, :, 0])
                
                if verbose:
                    print f + '_pert_comp' + row + ' loaded'
                    print 'Shape: ' + self.perts[f + row].shape
                    
                # raw complex result from calling 'read_field', which peaks at I-coil center
                # these need to be processed to get the correct toroidal angle and amplitudes
        
    
    # -----------------------------------------------------------------------
    
    def plot_perts(self, field, ntor, torloc_u = 30, torloc_l = 30, linfac = -2.5 * 4 / np.pi, norm = 1, clims = None,
                   add_eq = False, no_labels = False, plot_wall = True,
                   plot_Icoils = True, cmap = 'bwr', plot_psin = 0.94, pellet_ports = None,
                   exclude_outside = True, plot_gradient = False, save_fig_loc = None):
        """
        Plot the cross section of the perturbed quantities

        Inputs
          field       Selection of which perturbations to plot ('te', 'ti', 'p', 'pe' or 'ne')
          torloc      Toroidal angle for cross-section plot (in DIII-D toroidal coordinates)
          linfac      Linear factor to multiply the calculated perturbation by.
                      Scale up by 4/pi times I-coil current in kA
          norm        Value to normalize the qunatities by (i.e., 1e18 for density)
          clims       Should be given as a single number if add_eq is False, and the color map will be
                      symmetric about that
          cmap        Colormap choice ('bwr' good if not adding equilibrium, 'magma' if yes)
          add_eq      Adds the equilibrium values to the perturbed quantities
          ntor        Toroidal mode number of the applied perturbation
          pellet_ports List of pellet injection trajectories to plot ("HFS_Upper')
          plot_gradient Plots the radial gradient (meant for finding predicted islands)
          save_fig_loc  If not None, saves the figure using the given file location (as .png)
          torloc_*     Toroidal location (in degrees) of upper and lower I-coil rows (_u and _l)
                       in DIII-D coordinates
        """
        phi_in_l = np.deg2rad(torloc_l)
        phi_in_u = np.deg2rad(torloc_u)
        
        perts = np.zeros(self.perts[field + '_l'].shape)
        for row in ['_l', '_u']:
            perts_comp = self.perts[field + row] / norm  # complex expression of perturbations
            amp = linfac * np.sqrt(np.imag(perts_comp) ** 2 + np.real(perts_comp) ** 2)
            phase = np.arctan2(np.imag(perts_comp), np.real(perts_comp))
            
            # sum contributions of upper and lower rows
            if row == '_l':
                perts += amp * np.real(np.cos(phase + phi_in_l * ntor))
            else:
                perts += amp * np.real(np.cos(phase + phi_in_u * ntor))
        
        if exclude_outside:
            perts[self.psin2d > 1] = np.nan
            perts = np.ma.masked_invalid(perts)
        
        z_grid, r_grid = np.meshgrid(self.z, self.r)
        
        if add_eq:
            total = self.eqs[field] + perts
        else:
            total = perts.copy()
        
        f = plt.figure()
        if plot_gradient:
            # assume that dr and dz are constant
            dT_dr = np.zeros([len(self.r), len(self.z)])
            dT_dz = np.zeros([len(self.r), len(self.z)])
            
            dT_dr[1:-1, :] = (total[:-2, :] - total[2:, :]) / self.r[1] - self.r[0]
            dT_dz[:, 1:-1] = (total[:, :-2] - total[:, 2:]) / self.z[1] - self.z[0]
            dT_dr[0, :] = dT_dr[1, :]
            dT_dr[-1, :] = dT_dr[-2, :]
            dT_dz[:, 0] = dT_dz[:, 1]
            dT_dz[:, -1] = dT_dz[:, -2]
            
            grad = np.sqrt(dT_dr ** 2 + dT_dz ** 2)
            if exclude_outside:
                grad[self.psin2d > 1] = np.nan
                grad = np.ma.masked_invalid(grad)
            
            plt.pcolormesh(r_grid, z_grid, grad, cmap = cmap)
        else:
            plt.pcolormesh(r_grid, z_grid, total, cmap = cmap)
        
        cb = plt.colorbar()
        if clims is not None:
            if add_eq:
                plt.clim(clims)
            else:
                plt.clim([-clims, clims])
                # cb.set_ticks(range(-clims, clims + 1))
        
        self._set_plot_defaults(f, no_labels = no_labels, plot_wall = plot_wall,
                                plot_Icoils = plot_Icoils)
        if not no_labels:
            plt.title(field)
        
        if plot_psin is not None:
            psin = self.psin2d.copy()
            psin[self.psin2d > 1.01] = np.nan
            psin_masked = np.ma.masked_invalid(psin)
            cntr = plt.contour(r_grid, z_grid, psin_masked, [1], colors = '0')
            cntr = plt.contour(r_grid, z_grid, -psin_masked, [-plot_psin], colors = 'w')
            
            # plt.clabel(cntr, plot_psin, fontsize = 10, fmt = '%.12f', inline = 1)
        
        if pellet_ports is not None:
            for prt in pellet_ports:
                R_traj, Z_traj = get_pellet_traj(prt)
                plt.plot(R_traj, Z_traj, 'g', lw = 2)
                
        if save_fig_loc is not None:
            plt.savefig(save_fig_loc + '.png', format = 'png')
    
    # -----------------------------------------------------------------------
    
    def plot_eq(self, field, norm = 1, clims = None, no_labels = False, plot_wall = True, plot_Icoils = True,
                cmap = 'inferno'):
        """
        Plot equilibrium quantities in the cross section

          cmap    Set to the colormap that you want, or choose 'contour' for a contour plot
        """
        r_grid, z_grid = np.meshgrid(self.r, self.z)
        
        eq = self.eqs[field].copy() / norm
        
        f = plt.figure()
        
        if cmap == 'contour':
            eq[self.psin2d > 1.01] = np.nan
            eq = np.ma.masked_invalid(eq)
            
            cntr = plt.contour(r_grid, z_grid, eq, colors = '0')
            plt.clabel(cntr, fontsize = 10, fmt = '%i', inline = 1)
        else:
            eq[self.psin2d > 1] = np.nan
            eq = np.ma.masked_invalid(eq)
            plt.pcolormesh(r_grid, z_grid, eq, cmap = cmap)
            cb = plt.colorbar()
        
        if clims is not None:
            plt.clim(clims)
        
        self._set_plot_defaults(f, no_labels = no_labels, plot_wall = plot_wall,
                                plot_Icoils = plot_Icoils)

    # -----------------------------------------------------------------------

    def plot_xphase(self, fields = ['ne', 'te'], ul_phase = 0, linfac = -2.5 * 4 / np.pi, amp_mask = None,
                   clims = 180, no_labels = False, plot_wall = True,
                   plot_Icoils = True, cmap = None, exclude_outside = True, save_fig_loc = None, ret_vars = False):
        """
        Plot the cross phase of the perturbed quantities

        Inputs
          field       Selection of which perturbations to plot ('te', 'ti', 'p', 'pe' or 'ne')
          ul_phase    Difference between upper and lower coil phases (not ready for !=0 yet)
          linfac      Linear factor to multiply the calculated perturbation by.
                      Scale up by 4/pi times I-coil current in kA
          clims       Should be given as a single number, and the color map will be
                      symmetric about that
          cmap        Colormap choice (default is custom circular one, good for phase)
          save_fig_loc  If not None, saves the figure using the given file location (as .png)
        """
        
        if ul_phase != 0:
            print "code not ready for phase not = 0 yet"
            return
        
        perts_comp1 = self.perts[fields[0] + '_l'] + self.perts[fields[0] + '_u']
        amp = linfac * np.sqrt(np.imag(perts_comp1) ** 2 + np.real(perts_comp1) ** 2)
        phase1 = np.arctan2(np.imag(perts_comp1), np.real(perts_comp1))

        perts_comp2 = self.perts[fields[1] + '_l'] + self.perts[fields[1] + '_u']
        phase2 = np.arctan2(np.imag(perts_comp2), np.real(perts_comp2))
        
        xphase = phase2 - phase1
        
        xphase[xphase < -np.pi] += 2*np.pi
        xphase[xphase > np.pi] -= 2*np.pi

        # perts = np.zeros(self.perts[field + '_l'].shape, dtype = complex)
        #
        # perts_comp = self.perts[field + '_l'] + self.perts[field + '_u']
        #
        # amp = linfac * np.sqrt(np.imag(perts_comp) ** 2 + np.real(perts_comp) ** 2)
        # phase = np.arctan2(np.imag(perts_comp), np.real(perts_comp))
        #
        # for row in ['_l', '_u']:
        #     perts_comp = self.perts[field + row]  # complex expression of perturbations
        #     amp = linfac * np.sqrt(np.imag(perts_comp) ** 2 + np.real(perts_comp) ** 2)
        #     phase = np.arctan2(np.imag(perts_comp), np.real(perts_comp))
        #
        #     # sum contributions of upper and lower rows
        #     if row == '_l':
        #         perts += amp * np.real(np.cos(phase + phi_in_l * ntor))
        #     else:
        #         perts += amp * np.real(np.cos(phase + phi_in_u * ntor))
        
    
        if exclude_outside:
            xphase[self.psin2d > 1] = np.nan
            
            if amp_mask is None:
                xphase = np.ma.masked_invalid(xphase)
            
        if amp_mask is not None:
            xphase[np.abs(amp) < amp_mask] = np.nan
            xphase = np.ma.masked_invalid(xphase)
    
        phase_deg = np.rad2deg(xphase)
        z_grid, r_grid = np.meshgrid(self.z, self.r)
    
        f = plt.figure()
        if cmap is None:
            plt.pcolormesh(r_grid, z_grid, phase_deg, cmap = make_segmented_cmap(), vmin = -np.pi, vmax = np.pi)
        else:
            plt.pcolormesh(r_grid, z_grid, phase_deg, cmap = cmap)
    
        cb = plt.colorbar()
        if clims is not None:
            plt.clim([-clims, clims])
                # cb.set_ticks(range(-clims, clims + 1))

        if not no_labels:
            plt.title('$\phi$_' + fields[1] + ' - $\phi$_' + fields[0])
            
        self._set_plot_defaults(f, no_labels = no_labels, plot_wall = plot_wall,
                                plot_Icoils = plot_Icoils)
    
        if save_fig_loc is not None:
            plt.savefig(save_fig_loc + '.png', format = 'png')
            
        if ret_vars: return amp, phase1, phase2

    # -----------------------------------------------------------------------
    
    def plot_psin(self, clims = [0, 1], no_labels = False, plot_wall = True, plot_Icoils = True,
                  cmap = 'inferno_r', levels = [0.2, 0.4, 0.6, 0.8, 0.940, 1]):
        """
        Plot psin in the cross section

          cmap    Set to the colormap that you want, or choose 'contour' for a contour plot
        """
        z_grid, r_grid = np.meshgrid(self.z, self.r)
        
        psin = self.psin2d.copy()
        
        f = plt.figure()
        
        if cmap == 'contour':
            psin[self.psin2d > 1.01] = np.nan
            psin_masked = np.ma.masked_invalid(psin)
            
            if levels is None:
                levels = np.linspace(0, 1, 6)
            cntr = plt.contour(r_grid, z_grid, psin_masked, levels, colors = '0')
            plt.clabel(cntr, levels, fontsize = 10, fmt = '%.12f', inline = 1)
        else:
            psin[self.psin2d > 1] = np.nan
            psin_masked = np.ma.masked_invalid(psin)
            plt.pcolormesh(r_grid, z_grid, psin_masked, cmap = cmap)
            cb = plt.colorbar()
        
        if clims is not None:
            plt.clim(clims)
        
        self._set_plot_defaults(f, no_labels = no_labels, plot_wall = plot_wall,
                                plot_Icoils = plot_Icoils)

    # -----------------------------------------------------------------------
        
    def make_UL_phase_scan_gif(self, field = 'ne', nphases = 20, ntor=2, torloc = 30, linfac = -2.5 * 4 / np.pi,
                               norm = 1, clims = None, cmap = 'bwr', plot_psin = 0.94,
                               save_folder = None):
        """
        Make a series of figures using 'plot_perts'
        """
        
        ul_phases = np.linspace(0, 2*np.pi, nphases)
        
        
        for i, ph in enumerate(ul_phases):
            if save_folder is None:
                save_fig_loc = None
            else:
                if save_folder[-1] != '/':
                    save_folder += '/'
                    
                save_fig_loc = save_folder + '_' + str(i)
            
            
            self.plot_perts(field, ntor = ntor, torloc = torloc, phase_u_l = ph, linfac = linfac, norm = norm,
                            clims = clims, add_eq = False, no_labels = False, plot_wall = True,
                            plot_Icoils = True, cmap = cmap, plot_psin = plot_psin, pellet_ports = None,
                            exclude_outside = True, plot_gradient = False, save_fig_loc = save_fig_loc)


    # -----------------------------------------------------------------------
    def plot_xphase_1D(self, zloc = 0, r_range = [2.2, 2.276], ul_phase = 0, plot_eq = True,
                       linfac = -2.5 * 4 / np.pi):
        """
        Plot the cross phase of the perturbed quantities

        Inputs
          field       Selection of which perturbations to plot ('te', 'ti', 'p', 'pe' or 'ne')
          ul_phase    Difference between upper and lower coil phases (not ready for !=0 yet)
          linfac      Linear factor to multiply the calculated perturbation by.
                      Scale up by 4/pi times I-coil current in kA
        """
        
        fields = ['ne', 'te']  # Might as well hard code this, not looking at anything else and
                               #  that ingrained assumption might break things if we try to change it
    
        if ul_phase != 0:
            print "code not ready for phase not = 0 yet"
            return
        
        dr = np.min(np.abs(np.diff(self.r)))
        z_ind = np.argmin(np.abs(self.z - zloc))
        r_inds = np.where((self.r >= r_range[0] - dr) & (self.r <= r_range[1] + dr))
        Rlocs = self.r[r_inds]
    
        perts_comp1 = self.perts[fields[0] + '_l'][r_inds, z_ind][0,:] + \
                         self.perts[fields[0] + '_u'][r_inds, z_ind][0,:]
        # perts_comp1 = perts_comp1_2D[0,:]
        amp1 = linfac * np.sqrt(np.imag(perts_comp1) ** 2 + np.real(perts_comp1) ** 2)
        phase1 = np.arctan2(np.imag(perts_comp1), np.real(perts_comp1))
        eqs1 = self.eqs[fields[0]][r_inds, z_ind][0,:]
    
        perts_comp2 = self.perts[fields[1] + '_l'][r_inds, z_ind][0,:] + \
                         self.perts[fields[1] + '_u'][r_inds, z_ind][0,:]
        # perts_comp2 = perts_comp2_2D[0, :]
        amp2 = linfac * np.sqrt(np.imag(perts_comp2) ** 2 + np.real(perts_comp2) ** 2)
        phase2 = np.arctan2(np.imag(perts_comp2), np.real(perts_comp2))
        eqs2 = self.eqs[fields[1]][r_inds, z_ind][0,:]
    
        xphase = phase2 - phase1
    
        xphase[xphase < -np.pi] += 2 * np.pi
        xphase[xphase > np.pi] -= 2 * np.pi
        
        # perts = np.zeros(self.perts[field + '_l'].shape, dtype = complex)
        #
        # perts_comp = self.perts[field + '_l'] + self.perts[field + '_u']
        #
        # amp = linfac * np.sqrt(np.imag(perts_comp) ** 2 + np.real(perts_comp) ** 2)
        # phase = np.arctan2(np.imag(perts_comp), np.real(perts_comp))
        #
        # for row in ['_l', '_u']:
        #     perts_comp = self.perts[field + row]  # complex expression of perturbations
        #     amp = linfac * np.sqrt(np.imag(perts_comp) ** 2 + np.real(perts_comp) ** 2)
        #     phase = np.arctan2(np.imag(perts_comp), np.real(perts_comp))
        #
        #     # sum contributions of upper and lower rows
        #     if row == '_l':
        #         perts += amp * np.real(np.cos(phase + phi_in_l * ntor))
        #     else:
        #         perts += amp * np.real(np.cos(phase + phi_in_u * ntor))
    
        # Get the pressure gradient to compare with the cross phase
        
        pres = self.eqs['p'][r_inds, z_ind][0,:]
        grad_p = np.abs(np.diff(pres) / np.diff(Rlocs))
        Rlocs_mid = Rlocs[:-1] + np.diff(Rlocs) / 2
        
        # Find where psi_n = 1
        
        psin_1D = self.psin2d[r_inds, z_ind][0,:]
        psin_func = scinter.interp1d(psin_1D, Rlocs)
        Rloc_sep = psin_func(1)
    
        # Plot everything

        f, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex = True, figsize=(8,10))
        
        ax0.plot(Rlocs, np.rad2deg(xphase), 'k', lw = 3, zorder = 5)
        ax0.set_title('$\phi$_' + fields[0] + ' - ' + '$\phi$_' + fields[1] +
                      ' (Z = {0:.3f}'.format(self.z[z_ind]) + ' m)')
        ax0.set_ylim([-180, 180])
        ax0.set_yticks(np.arange(-180, 181, 90))
        ax0.set_yticklabels(['-180$^\circ$', '-90$^\circ$', '0$^\circ$', '90$^\circ$', '180$^\circ$'])
        
        ax1.plot(Rlocs_mid, grad_p/1e3, 'g', lw=2, zorder = 5)
        # ax1.set_ylabel('|$\\nabla$P| (kPa / m)')
        ax1.set_ylabel('|$\\nabla$P|')
        ax1.set_yticks(range(0, 650, 200))
        
        ax2.semilogy(Rlocs, np.abs(amp1), 'b', lw=2, label = 'Pert n = 3', zorder = 5)
        # ax2.set_ylabel('$\delta' + fields[0])
        ax2.set_ylabel('$\delta$n$_e$')

        ax3.semilogy(Rlocs, np.abs(amp2), 'r', lw=2, label = 'Pert n = 3', zorder = 5)
        # ax3.set_ylabel('$\delta$' + fields[1])
        ax3.set_ylabel('$\delta$T$_e$')
        ax3.set_xlabel('R (m)')
        
        for i, a in enumerate([ax2, ax3]):
            if fields[i] == 'ne':
                a.set_ylim([2e16, 4e19])
            elif fields[i] == 'te':
                a.set_ylim([5, 3e3])
            
        for a in [ax0, ax1, ax2, ax3]:
            a.grid('on')
            y_lims = a.get_ylim()
            a.plot([Rloc_sep, Rloc_sep], y_lims, '--k', lw = 3, zorder = 1)
            
        if plot_eq:
            ax2.semilogy(Rlocs, eqs1, '--b', lw=1, label = 'n = 0', zorder = 4)
            ax2.legend(loc = 'center left', fontsize = 14)
            ax3.semilogy(Rlocs, eqs2, '--r', lw=1, label = 'n = 0', zorder = 4)
            ax3.legend(loc = 'center left', fontsize=14)
            
        ax3.set_xlim(r_range)
            
        plt.tight_layout(pad = 0.5)
        plt.show(block = False)
    
            
    
    # -----------------------------------------------------------------------
    
    def _set_plot_defaults(self, f, no_labels = False, plot_wall = True, plot_Icoils = True):
        
        if plot_wall:
            R_wall, Z_wall = get_wall()
            plt.plot(R_wall, Z_wall, '-k', zorder = 1)
        
        if plot_Icoils:
            R_IL, Z_IL, R_IU, Z_IU = get_Icoils()
            plt.plot(R_IL, Z_IL, '-k', lw = 4, zorder = 2)
            plt.plot(R_IU, Z_IU, '-k', lw = 4, zorder = 2)
        
        plt.axis('equal')
        f.axes[0].set_xlim(min(self.r) - 0.1, max(self.r) + 0.1)
        f.axes[0].set_ylim(min(self.z) - 0.12, max(self.z) + 0.04)
        if no_labels:
            plt.xticks([1, 1.5, 2, 2.5], ['', '', '', ''])
            plt.yticks([-1, -0.5, 0, 0.5, 1], ['', '', '', '', ''])
        else:
            plt.xticks([1, 1.5, 2, 2.5])
            plt.yticks([-1, -0.5, 0, 0.5, 1])
            plt.xlabel('R (m)')
            plt.ylabel('Z (m)')
        
        plt.tight_layout()
        plt.show(block = False)

#-----------------------------------------------------------------------


class M3DC1_OMP(object):
    def __init__(self, sav_fileloc = '/home/wilcoxr/m3dc1/profile_savs/170007/torrad_170007_3625.sav'):
        """
        Requires previously saved IDL .sav file containing profile variables
        I generate this with 'load_m3dc1.pro', here on venus: /u/wilcoxr/idl/load_m3dc1.pro
           or 'load_m3dc1_single_n.pro'
        """
        
        idlsave_dict = si.readsav(sav_fileloc, python_dict = True, verbose = False)
        
        # fields = ['te','den','ti','pe','p'] (usually ordered like this in arrays)
        
        
        self.fields = idlsave_dict['fields']  # fields that are in the IDL sav file, and which index corresponds to which
        self.r = idlsave_dict['r']  # grid locations for original 'read_field' call
        self.z = idlsave_dict['z']
        self.rv = idlsave_dict['rv']  # grid locations at outboard midplane
        self.zv = idlsave_dict['zv']
        self.psinrv = idlsave_dict['psinrv']  # psi at outboard midplane
        self.eqs = idlsave_dict['eqs']  # equilibrium fields (field, surf) (5, 701)
        self.perts_tot = idlsave_dict['perts_tot']  # perterbed quantities at BES, toroidal array
                                                    #  (tor angle, field, surf) (360, 5, 701)
        # self.pertsl_tot = idlsave_dict['pertsl_tot']  # perterbed gradient scale lengths (360, 5, 701) (maybe add this)
        
        self.phi = idlsave_dict['phi_machine_tot']
        
        self.dens_ind = np.where(self.fields == 'den')[0][0]
        self.te_ind = np.where(self.fields == 'te')[0][0]
        self.ti_ind = np.where(self.fields == 'ti')[0][0]
        
        # nprof = len(self.perts_tot[0, 0, :])
    
    # -----------------------------------------------------------------------
    
    def plot_dens_pert_OMP(self, psi_range = [0.924, 1.009], fig_size = (10, 6),
                           plot_diagnostics = True, plot_labels = True, lw_border = 2):
        """
        Plot the surface (Te) and density contours (NOT color contours) at the outboard midplane
        as a function of radial and toroidal location

        Inputs:
          plot_diagnostics
        """
        # npsi = len(self.psinrv)
        phi = self.phi
        
        rad_grid, tor_grid = np.meshgrid(self.psinrv, phi)  # toroidal locs X surfaces
        
        dens_pert = self.perts_tot[:, self.dens_ind, :]
        
        fig = plt.figure(figsize = fig_size)
        plt.pcolormesh(tor_grid, rad_grid, dens_pert, cmap = 'RdBu')
        ax = fig.axes[0]
        
        if plot_diagnostics:
            plt.text(145, 1, '$\\downarrow$', fontsize = 30, color = 'k', fontweight = 'bold',
                     horizontalalignment = 'center', verticalalignment = 'bottom')
            plt.text(255, 1, '$\\downarrow$', fontsize = 30, color = 'k', fontweight = 'bold',
                     horizontalalignment = 'center', verticalalignment = 'bottom')
            if plot_labels:
                plt.text(152, 1, 'BES', fontsize = 20, color = 'k', fontweight = 'bold',
                         horizontalalignment = 'left', verticalalignment = 'bottom')
                plt.text(262, 1, 'Reflectometer', fontsize = 20, color = 'k', fontweight = 'bold',
                         horizontalalignment = 'left', verticalalignment = 'bottom')
        
        plt.ylim(psi_range)
        plt.xlim([0, 359])
        
        if plot_labels:
            plt.xlabel('$\\phi$ (deg)', fontsize = 20)
            plt.ylabel('$\\psi_N$', fontsize = 20)
        plt.xticks(np.arange(0, 360, 60))
        plt.tight_layout()
        
        # ----------------------------------------
        # Set major and minor ticks to look nice
        
        ymajorLocator = MultipleLocator(0.05)
        yminorLocator = MultipleLocator(0.01)
        
        xmajorLocator = MultipleLocator(60)
        xminorLocator = MultipleLocator(10)
        
        ax.xaxis.set_major_locator(xmajorLocator)
        ax.xaxis.set_minor_locator(xminorLocator)
        
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_minor_locator(yminorLocator)
        
        plt.tick_params(which = 'both', width = lw_border)
        plt.tick_params(which = 'major', length = 7)
        plt.tick_params(which = 'minor', length = 4)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(lw_border)
            
        plt.show(block = False)
    
    # -----------------------------------------------------------------------
    
    def plot_surf_contours_OMP(self, psi_range = [0.93, 1.001], dpsi_label = 0.01, dpsi_cont = 0.01, dR = 0.005,
                               ndens = 30, nTe = 30, fig_size = (10, 6), lw_border = 2,
                               shift_60 = False, clim_max = 120, R_range = [2.25, 2.276],
                               plot_eq = False, plot_Te = False, plot_aLTe = False, plot_psi_contours = False,
                               R_not_psi = False, plot_diagnostics = False, plot_labels = True, plot_Ti = False,
                               cmap_d = 'inferno', cmap_t = 'magma'):
        """
        Plot the density perturbations at the outboard midplane
        as a function of radial and toroidal location

        Radial positions are re-mapped onto constant Te-contours by having R_not_psi = False

        Inputs:
          plot_diagnostics
          shift_60           **doesn't work yet**
          R_not_psi          Plot y-axis variable as major radius rather than flux
                             Leaving this as 'False' maps contours onto constant-Te
        """
        Te_thresh = 1  # Don't plot anything below this
        
        npsi = len(self.psinrv)
        ntor = len(self.perts_tot[:, 0, 0])
        
        phi = np.linspace(0, 360, ntor + 1)[:-1]
        if shift_60:
            perts_tot = np.roll(self.perts_tot, 60, axis = 0)
        else:
            perts_tot = self.perts_tot
        
        psi_grid_orig, tor_grid = np.meshgrid(self.psinrv, phi)  # toroidal locs X surfaces
        R_grid, tor_grid = np.meshgrid(self.rv, phi)  # toroidal locs X surfaces
        dens_sym, _ = np.meshgrid(self.eqs[self.dens_ind, :], phi)  # toroidal locs X surfaces
        Te_sym, _ = np.meshgrid(self.eqs[self.te_ind, :], phi)  # toroidal locs X surfaces
        
        dens = dens_sym + perts_tot[:, self.dens_ind, :]  # total density (equilib + perturbed)
        Te = Te_sym + perts_tot[:, self.te_ind, :]  # total Te (equilib + perturbed)
        Te[np.where(Te < Te_thresh)] = Te_thresh  # don't plot the negative Te contours in the SOL
        
        R_grid_mid = R_grid[:, :-1] + np.diff(R_grid, axis = 1) / 2
        
        # Map grid points to psi values using constant-Te contours
        fitted_Te_inds = np.where(self.psinrv < 1.1)[0]
        f_Te_to_psi = scinter.interp1d(self.eqs[self.te_ind, fitted_Te_inds], self.psinrv[fitted_Te_inds],
                                       fill_value = 'extrapolate')
        # fill_value = [1, np.max(self.eqs[self.te_ind,:])])
        psi_grid = f_Te_to_psi(Te)
        self.psi_grid = psi_grid
        psi_conts = np.arange(1, psi_range[0] - dpsi_cont * 1.1, -dpsi_cont)[1:]  # psi contours, add [::-1] ascending
        
        if R_not_psi:
            y_grid = R_grid
            y_grid_orig = R_grid
            y_grid_mid = R_grid_mid
            y_grid_orig_mid = R_grid_mid
            f_psi_to_R_eq = scinter.interp1d(self.psinrv[fitted_Te_inds], self.rv[fitted_Te_inds])
            y_LCFS = f_psi_to_R_eq(1.0)
            y_range = R_range
            # y_range = f_psi_to_R_eq(psi_range)
            
            # Find LCFS as a function of toroidal angle, to plot contour
            if plot_eq:
                R_LCFS_3D = np.ones(ntor) * f_psi_to_R_eq(1)
            else:
                R_LCFS_3D = np.zeros(ntor)
                # R_dash_3D = np.zeros([ntor, len(psi_conts)])
                for i, phi_tmp in enumerate(phi):
                    f_psi_to_R = scinter.interp1d(psi_grid[i, fitted_Te_inds], self.rv[fitted_Te_inds])
                    R_LCFS_3D[i] = f_psi_to_R(1)
                    
                    # if plot_psi_contours:
                    #     for p_ind, p in enumerate(psi_conts):
                    #         R_dash_3D[i, p_ind] = f_psi_to_R(p)
        else:
            y_grid = psi_grid
            y_grid_orig = psi_grid_orig
            psi_grid_mid = psi_grid[:, :-1] + np.diff(psi_grid, axis = 1) / 2
            y_grid_orig_mid = psi_grid_orig[:, :-1] + np.diff(psi_grid_orig, axis = 1) / 2
            y_grid_mid = psi_grid_mid
            y_LCFS = 1
            y_range = psi_range
        
        # Calculate gradient scale lengths
        dens_mid = dens[:, :-1] + np.diff(dens, axis = 1) / 2.
        aLne_mid = a_minor * np.abs((np.diff(dens, axis = 1) / np.diff(self.rv))) / dens_mid
        # aLne = np.append(aLne_mid, np.array([aLne_mid[:, -1]]).T, axis = 1)
        
        
        # grad_n_mid = np.diff(dens, axis = 1) / np.diff(self.r_locs)  # on the midpoint grid
        # grad_Te_mid = np.diff(Te, axis = 1) / np.diff(self.r_locs)  # on the midpoint grid
        #
        # # interpolate at each toroidal location
        # # (probably could do this with interp1d and axis=1, but it didn't work right away)
        # grad_n_rlocs = np.zeros([ntor, npsi - 2])
        # grad_Te_rlocs = np.zeros([ntor, npsi - 2])
        # for i, tor in enumerate(phi):
        #     f_gradn = scinter.interp1d(r_locs_mid, grad_n_mid[i, :], kind = 'cubic')
        #     grad_n_rlocs[i, :] = f_gradn(self.r_locs[1:-1]) # grad_n at the r_loc locations
        #
        #     if plot_aLTe:
        #         f_gradTe = scinter.interp1d(r_locs_mid, grad_Te_mid[i, :], kind = 'cubic')
        #         grad_Te_rlocs[i, :] = f_gradTe(self.r_locs[1:-1]) # grad_Te at the r_loc locations
        #
        # # try doing without loop (doesn't work)
        # # f_gradn = scinter.interp1d(r_locs_mid, grad_n_mid, kind = 'cubic', axis = 1)
        # # grad_n_rlocs = f_gradn(self.r_locs[1:-1]) # grad_n at the r_loc locations
        #
        # grad_n_almost = np.append(np.array([grad_n_mid[:, 0]]).T, grad_n_rlocs, axis = 1)  # copy endpoints
        # grad_n = np.append(grad_n_almost, np.array([grad_n_mid[:, -1]]).T, axis = 1)  # copy endpoints
        # aLne = a_minor * grad_n / dens
        
        
        # Plot density and density gradient
        
        fig = plt.figure(figsize = fig_size)
        if plot_eq:
            dens_sym_mid = dens_sym[:, :-1] + np.diff(dens_sym, axis = 1) / 2.
            aLne_eq_mid = a_minor * np.abs((np.diff(dens_sym, axis = 1) / np.diff(self.rv))) / dens_sym_mid
            # aLne_eq = np.append(aLne_eq_mid, np.array([aLne_eq_mid[:, -1]]).T, axis = 1)
            
            plt.pcolormesh(tor_grid[:, :-1], y_grid_orig_mid, aLne_eq_mid, cmap = cmap_d)
        
        else:
            plt.pcolormesh(tor_grid[:, :-1], y_grid_mid, aLne_mid, cmap = cmap_d)
        
        plt.colorbar()
        if clim_max is not None:
            plt.clim([0, clim_max])
        
        # Plot contours of ne (and Te if wanted)
        if plot_Te:
            if plot_eq:
                ct_Te = plt.contour(tor_grid, y_grid_orig, -Te_sym, nTe, colors = 'w')  # negative makes it dashed
            else:
                ct_Te = plt.contour(tor_grid, y_grid, -Te, nTe, colors = 'w')  # negative makes it dashed
        
        if plot_psi_contours:
            if R_not_psi:
                if plot_eq:
                    # for p in psi_conts:
                    #     R_dash = f_psi_to_R_eq(p)
                    #     plt.plot([0,359], [R_dash, R_dash], '--w', lw = 1)
                    ct_psi = plt.contour(tor_grid, y_grid, -psi_grid_orig, colors = 'w', levels = -psi_conts)
                else:
                    # for p_ind, p in enumerate(psi_conts):
                    #     plt.plot(phi, R_dash_3D[:, p_ind], '--w', lw=1)
                    ct_psi = plt.contour(tor_grid, y_grid, -psi_grid, colors = 'w', levels = -psi_conts)
                    
                    # if dpsi_cont < 0.01:
                    #     plt.clabel(ct_psi, fontsize = 12, fmt='%1.3f', inline = 0) #, color = 'w'
                    # else:
                    #     plt.clabel(ct_psi, fontsize = 12, fmt='%1.2f', inline = 0) #, color = 'w'
            
            else:
                for p in psi_conts:
                    plt.plot([0, 359], [p, p], '--w', lw = 1)  # No need to label, already on y-axis
        
        if plot_eq:
            ct_dens = plt.contour(tor_grid, y_grid_orig, dens_sym, ndens, colors = 'r')
        else:
            ct_dens = plt.contour(tor_grid, y_grid, dens, ndens, colors = 'r')
        
        if plot_diagnostics:
            plt.plot([145, 145], y_range, '--c', lw = 3)
            plt.plot([255, 255], y_range, '--c', lw = 3)
            if plot_labels:
                text_y_loc = y_range[0] + (y_range[1] - y_range[0]) / y_range[1] * 0.2
                plt.text(145, text_y_loc, 'BES', fontsize = 20, color = 'c', fontweight = 'bold',
                         horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
                plt.text(255, text_y_loc, 'Reflectometer', fontsize = 20, color = 'c', fontweight = 'bold',
                         horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
        
        if R_not_psi:
            # plt.yticks(np.arange(np.around(y_LCFS, decimals = 2), y_range[0] - dR/2, -dR))
            plt.yticks(np.arange(np.ceil(y_range[0] / dR) * dR, np.floor(y_range[1] / dR) * dR + dR / 2, dR))
            if plot_psi_contours or plot_Te:
                plt.plot(phi, R_LCFS_3D, '--w', lw = 3)
        else:
            plt.plot([0, 359], [1, 1], '--w', lw = 3)
            plt.yticks(np.arange(1., psi_range[0] - dpsi_label / 2, -dpsi_label))
        
        plt.ylim(y_range)
        plt.xlim([0, 359])
        
        if plot_labels:
            plt.xlabel('$\\phi$ (deg)', fontsize = 20)
            if R_not_psi:
                plt.ylabel('R (m)', fontsize = 20)
            else:
                plt.ylabel('$\\psi_N$', fontsize = 20)
        
        plt.xticks(np.arange(0, 360, 60))
        plt.tight_layout()
        
        plt.tick_params(which = 'both', width = lw_border)
        plt.tick_params(which = 'major', length = 7)
        plt.tick_params(which = 'minor', length = 4)
        for axis in ['top', 'bottom', 'left', 'right']:
            fig.axes[0].spines[axis].set_linewidth(lw_border)
        
        # Now plot a/LTe for comparison
        
        if plot_aLTe:
            
            fig = plt.figure(figsize = fig_size)
            if plot_eq:
                Te_sym_mid = Te_sym[:, :-1] + np.diff(Te_sym, axis = 1) / 2.
                aLTe_eq_mid = a_minor * np.abs((np.diff(Te_sym, axis = 1) / np.diff(self.rv))) / Te_sym_mid
                # aLTe_eq = np.append(aLTe_eq_mid, np.array([aLTe_eq_mid[:, -1]]).T, axis = 1)
                plt.pcolormesh(tor_grid[:, :-1], y_grid_orig_mid, aLTe_eq_mid, cmap = cmap_t)
            else:
                Te_mid = Te[:, :-1] + np.diff(Te, axis = 1) / 2.
                aLTe_mid = a_minor * np.abs((np.diff(Te, axis = 1) / np.diff(self.rv))) / Te_mid
                # aLTe = np.append(aLTe_short, np.array([aLTe_short[:, -1]]).T, axis = 1)
                plt.pcolormesh(tor_grid[:, :-1], y_grid_mid, aLTe_mid, cmap = cmap_t)
            plt.colorbar()
            if clim_max is not None:
                plt.clim([0, clim_max])
            
            if plot_Te:
                if plot_eq:
                    ct_Te = plt.contour(tor_grid, y_grid_orig, -Te_sym, nTe, colors = 'w')  # negative makes it dashed
                else:
                    ct_Te = plt.contour(tor_grid, y_grid, -Te, nTe, colors = 'w')  # negative makes it dashed
            
            if R_not_psi:
                plt.yticks(np.arange(np.around(y_LCFS, decimals = 2), y_range[0] - dR / 2, -dR))
                plt.plot(phi, R_LCFS_3D, '--w', lw = 3)
            else:
                plt.plot([0, 359], [1, 1], '--w', lw = 3)
                plt.yticks(np.arange(1., psi_range[0] - dpsi_label / 2, -dpsi_label))
            
            plt.ylim(y_range)
            plt.xlim([0, 359])
            
            if plot_diagnostics:
                plt.plot([145, 145], y_range, '--w', lw = 2)
                plt.plot([255, 255], y_range, '--w', lw = 2)
                if plot_labels:
                    text_y_loc = y_range[0] + (y_range[1] - y_range[0]) / y_range[1] * 0.2
                    plt.text(145, text_y_loc, 'BES', fontsize = 20, color = 'w', fontweight = 'bold',
                             horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
                    plt.text(255, text_y_loc, 'Reflectometer', fontsize = 20, color = 'w', fontweight = 'bold',
                             horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
            if plot_labels:
                plt.xlabel('$\\phi$ (deg)', fontsize = 20)
                plt.ylabel('$\\psi_N$', fontsize = 20)
            plt.xticks(np.arange(0, 360, 60))
            plt.tight_layout()
        
        if plot_Ti:
            
            fig = plt.figure(figsize = fig_size)
            
            Ti_sym, _ = np.meshgrid(self.eqs[self.ti_ind, :], phi)  # toroidal locs X surfaces
            
            Ti = Ti_sym + perts_tot[:, self.ti_ind, :]  # total Ti (equilib + perturbed)
            Ti[np.where(Ti < Te_thresh)] = Te_thresh  # don't plot the negative Ti contours in the SOL
            
            if plot_eq:
                plt.pcolormesh(tor_grid, y_grid_orig, Ti_sym, cmap = cmap_t)
            else:
                plt.pcolormesh(tor_grid, y_grid, Ti, cmap = cmap_t)
            
            plt.colorbar()
            if clim_max is not None:
                plt.clim([0, clim_max])
            
            if plot_Te:
                if plot_eq:
                    ct_Te = plt.contour(tor_grid, y_grid_orig, -Te_sym, nTe, colors = 'w')  # negative makes it dashed
                else:
                    ct_Te = plt.contour(tor_grid, y_grid, -Te, nTe, colors = 'w')  # negative makes it dashed
            
            # Overplot contours
            if plot_eq:
                ct_dens = plt.contour(tor_grid, y_grid_orig, Ti_sym, ndens, colors = 'r')
            else:
                ct_dens = plt.contour(tor_grid, y_grid, Ti, ndens, colors = 'r')
            
            if R_not_psi:
                plt.yticks(np.arange(np.around(y_LCFS, decimals = 2), y_range[0] - dR / 2, -dR))
                plt.plot(phi, R_LCFS_3D, '--w', lw = 3)
            else:
                plt.plot([0, 359], [1, 1], '--w', lw = 3)
                plt.yticks(np.arange(1., psi_range[0] - dpsi_label / 2, -dpsi_label))
            
            plt.ylim(y_range)
            plt.xlim([0, 359])
            
            if plot_diagnostics:
                plt.plot([145, 145], y_range, '--w', lw = 2)
                plt.plot([255, 255], y_range, '--w', lw = 2)
                if plot_labels:
                    text_y_loc = y_range[0] + (y_range[1] - y_range[0]) / y_range[1] * 0.2
                    plt.text(145, text_y_loc, 'BES', fontsize = 20, color = 'w', fontweight = 'bold',
                             horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
                    plt.text(255, text_y_loc, 'Reflectometer', fontsize = 20, color = 'w', fontweight = 'bold',
                             horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
            if plot_labels:
                plt.xlabel('$\\phi$ (deg)', fontsize = 20)
                plt.ylabel('$\\psi_N$', fontsize = 20)
            plt.xticks(np.arange(0, 360, 60))
            plt.tight_layout()

        plt.show(block = False)
            

#-----------------------------------------------------------------------


class M3DC1_OMP_157306(object):
    def __init__(self, sav_fileloc = '/home/wilcoxr/m3dc1/profile_savs/157306/157306_m3dc1_n3_nf_allfields.sav'):
        """
        Requires previously saved IDL .sav file containing profile variables
        I generate this with 'load_m3dc1.pro', here on venus: /u/wilcoxr/idl/load_m3dc1.pro
           or 'load_m3dc1_single_n.pro'
        """
        
        idlsave_dict = si.readsav(sav_fileloc, python_dict = True, verbose = False)
        
        # fields = ['te','den','ti','pe','p'] (usually ordered like this in arrays)
        

        self.fields = idlsave_dict['fields']  # fields that are in the IDL sav file, and which index corresponds to which
        self.r = idlsave_dict['r']  # grid locations (401)
        self.z = idlsave_dict['z']
        self.psinrv = idlsave_dict['psinrv']  # psi at outboard midplane? (701)
        self.eqs = idlsave_dict['eqs']  # equilibrium fields (field, surf) (5, 701)
        self.perts = idlsave_dict['perts']  # perterbed quantities at BES (0/60, field, surf) (2, 5, 701)
        self.perts_tot = idlsave_dict['perts_tot']  # perterbed quantities at BES, toroidal array
                                                    #  (tor angle, field, surf) (360, 5, 701)
        self.pertst = idlsave_dict['pertst'] # perterbed quantities at Thomson
        self.pertsr = idlsave_dict['pertsr'] # perterbed quantities at reflectometer
        self.eqst = idlsave_dict['eqst']
        self.pertsl_tot = idlsave_dict['pertsl_tot']  # perterbed gradient scale lengths (360, 5, 701)

        ntor = len(self.perts_tot[:, 0, 0])
        self.phi = np.linspace(0, 360, ntor + 1)[:-1]

        self.dens_ind = np.where(self.fields == 'den')[0][0]
        self.te_ind = np.where(self.fields == 'te')[0][0]
        self.ti_ind = np.where(self.fields == 'ti')[0][0]
        
        nprof = len(self.perts_tot[0, 0, :])
        # this is how it was defined in the code that generated the IDL sav files
        r_locs = np.linspace(0, nprof, nprof+1)[:-1] / (nprof-1)*0.3 + 2.1
        self.r_locs = np.array(r_locs)
        

    # -----------------------------------------------------------------------

    def plot_dens_pert_OMP(self, psi_range = [0.924, 1.009], fig_size = (10, 6),
                               plot_diagnostics = True, plot_labels = True, lw_border = 2):
        """
        Plot the surface (Te) and density contours (NOT color contours) at the outboard midplane
        as a function of radial and toroidal location

        Inputs:
          plot_diagnostics
        """
        # npsi = len(self.psinrv)
        phi = self.phi


        rad_grid, tor_grid = np.meshgrid(self.psinrv, phi)  # toroidal locs X surfaces

        dens_pert = self.perts_tot[:, self.dens_ind, :]

        fig = plt.figure(figsize = fig_size)
        plt.pcolormesh(tor_grid, rad_grid, dens_pert, cmap = 'RdBu')
        ax = fig.axes[0]

        if plot_diagnostics:
            plt.text(145, 1, '$\\downarrow$', fontsize = 30, color = 'k', fontweight = 'bold',
                     horizontalalignment = 'center', verticalalignment = 'bottom')
            plt.text(255, 1, '$\\downarrow$', fontsize = 30, color = 'k', fontweight = 'bold',
                     horizontalalignment = 'center', verticalalignment = 'bottom')
            if plot_labels:
                plt.text(152, 1, 'BES', fontsize = 20, color = 'k', fontweight = 'bold',
                         horizontalalignment = 'left', verticalalignment = 'bottom')
                plt.text(262, 1, 'Reflectometer', fontsize = 20, color = 'k', fontweight = 'bold',
                         horizontalalignment = 'left', verticalalignment = 'bottom')

        plt.ylim(psi_range)
        plt.xlim([0, 359])

        if plot_labels:
            plt.xlabel('$\\phi$ (deg)', fontsize = 20)
            plt.ylabel('$\\psi_N$', fontsize = 20)
        plt.xticks(np.arange(0, 360, 60))
        plt.tight_layout()

        # ----------------------------------------
        # Set major and minor ticks to look nice

        ymajorLocator = MultipleLocator(0.05)
        yminorLocator = MultipleLocator(0.01)
        
        xmajorLocator = MultipleLocator(60)
        xminorLocator = MultipleLocator(10)

        ax.xaxis.set_major_locator(xmajorLocator)
        ax.xaxis.set_minor_locator(xminorLocator)

        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_minor_locator(yminorLocator)

        plt.tick_params(which = 'both', width = lw_border)
        plt.tick_params(which = 'major', length = 7)
        plt.tick_params(which = 'minor', length = 4)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(lw_border)

    # -----------------------------------------------------------------------
    
    def plot_surf_contours_OMP(self, psi_range = [0.93, 1.001], dpsi_label = 0.01, dpsi_cont = 0.01, dR = 0.005,
                               ndens = 30, nTe = 30, fig_size=(10, 6), lw_border = 2,
                               shift_60 = True, clim_max = 120, R_range = [2.25, 2.276],
                               plot_eq = False, plot_Te = False, plot_aLTe = False, plot_psi_contours = False,
                               R_not_psi = False, plot_diagnostics = False, plot_labels = True, plot_Ti = False,
                               cmap_d = 'inferno', cmap_t = 'magma'):
        """
        Plot the density perturbations at the outboard midplane
        as a function of radial and toroidal location
        
        Radial positions are re-mapped onto constant Te-contours by having R_not_psi = False
        
        Inputs:
          plot_diagnostics
          shift_60           **doesn't work yet**
          R_not_psi          Plot y-axis variable as major radius rather than flux
                             Leaving this as 'False' maps contours onto constant-Te
        """
        Te_thresh = 1  # Don't plot anything below this
        
        npsi = len(self.psinrv)
        ntor = len(self.perts_tot[:, 0, 0])
        
        phi = np.linspace(0, 360, ntor + 1)[:-1]
        if shift_60:
            perts_tot = np.roll(self.perts_tot, 60, axis = 0)
        else:
            perts_tot = self.perts_tot
            

        psi_grid_orig, tor_grid = np.meshgrid(self.psinrv, phi)  # toroidal locs X surfaces
        R_grid, tor_grid = np.meshgrid(self.r_locs, phi)  # toroidal locs X surfaces
        dens_sym, _ = np.meshgrid(self.eqs[self.dens_ind, :], phi)  # toroidal locs X surfaces
        Te_sym, _ = np.meshgrid(self.eqs[self.te_ind, :], phi)  # toroidal locs X surfaces

        dens = dens_sym + perts_tot[:, self.dens_ind, :]  # total density (equilib + perturbed)
        Te = Te_sym + perts_tot[:, self.te_ind, :]  # total Te (equilib + perturbed)
        Te[np.where(Te < Te_thresh)] = Te_thresh  # don't plot the negative Te contours in the SOL

        R_grid_mid = R_grid[:, :-1] + np.diff(R_grid, axis=1) / 2
        
        # Map grid points to psi values using constant-Te contours
        fitted_Te_inds = np.where(self.psinrv < 1.1)[0]
        f_Te_to_psi = scinter.interp1d(self.eqs[self.te_ind, fitted_Te_inds], self.psinrv[fitted_Te_inds],
                                       fill_value = 'extrapolate')
                                       # fill_value = [1, np.max(self.eqs[self.te_ind,:])])
        psi_grid = f_Te_to_psi(Te)
        self.psi_grid = psi_grid
        psi_conts = np.arange(1, psi_range[0] - dpsi_cont * 1.1, -dpsi_cont)[1:]  # psi contours, add [::-1] ascending
        
        if R_not_psi:
            y_grid = R_grid
            y_grid_orig = R_grid
            y_grid_mid = R_grid_mid
            y_grid_orig_mid = R_grid_mid
            f_psi_to_R_eq = scinter.interp1d(self.psinrv[fitted_Te_inds], self.r_locs[fitted_Te_inds])
            y_LCFS = f_psi_to_R_eq(1.0)
            y_range = R_range
            # y_range = f_psi_to_R_eq(psi_range)
            
            # Find LCFS as a function of toroidal angle, to plot contour
            if plot_eq:
                R_LCFS_3D = np.ones(ntor) * f_psi_to_R_eq(1)
            else:
                R_LCFS_3D = np.zeros(ntor)
                # R_dash_3D = np.zeros([ntor, len(psi_conts)])
                for i, phi_tmp in enumerate(phi):
                    f_psi_to_R = scinter.interp1d(psi_grid[i, fitted_Te_inds], self.r_locs[fitted_Te_inds])
                    R_LCFS_3D[i] = f_psi_to_R(1)
                    
                    # if plot_psi_contours:
                    #     for p_ind, p in enumerate(psi_conts):
                    #         R_dash_3D[i, p_ind] = f_psi_to_R(p)
        else:
            y_grid = psi_grid
            y_grid_orig = psi_grid_orig
            psi_grid_mid = psi_grid[:, :-1] + np.diff(psi_grid, axis = 1) / 2
            y_grid_orig_mid = psi_grid_orig[:, :-1] + np.diff(psi_grid_orig, axis = 1) / 2
            y_grid_mid = psi_grid_mid
            y_LCFS = 1
            y_range = psi_range

            
        # Calculate gradient scale lengths
        dens_mid = dens[:, :-1] + np.diff(dens, axis = 1) / 2.
        aLne_mid = a_minor * np.abs((np.diff(dens, axis = 1) / np.diff(self.r_locs))) / dens_mid
        # aLne = np.append(aLne_mid, np.array([aLne_mid[:, -1]]).T, axis = 1)
        
        
        # grad_n_mid = np.diff(dens, axis = 1) / np.diff(self.r_locs)  # on the midpoint grid
        # grad_Te_mid = np.diff(Te, axis = 1) / np.diff(self.r_locs)  # on the midpoint grid
        #
        # # interpolate at each toroidal location
        # # (probably could do this with interp1d and axis=1, but it didn't work right away)
        # grad_n_rlocs = np.zeros([ntor, npsi - 2])
        # grad_Te_rlocs = np.zeros([ntor, npsi - 2])
        # for i, tor in enumerate(phi):
        #     f_gradn = scinter.interp1d(r_locs_mid, grad_n_mid[i, :], kind = 'cubic')
        #     grad_n_rlocs[i, :] = f_gradn(self.r_locs[1:-1]) # grad_n at the r_loc locations
        #
        #     if plot_aLTe:
        #         f_gradTe = scinter.interp1d(r_locs_mid, grad_Te_mid[i, :], kind = 'cubic')
        #         grad_Te_rlocs[i, :] = f_gradTe(self.r_locs[1:-1]) # grad_Te at the r_loc locations
        #
        # # try doing without loop (doesn't work)
        # # f_gradn = scinter.interp1d(r_locs_mid, grad_n_mid, kind = 'cubic', axis = 1)
        # # grad_n_rlocs = f_gradn(self.r_locs[1:-1]) # grad_n at the r_loc locations
        #
        # grad_n_almost = np.append(np.array([grad_n_mid[:, 0]]).T, grad_n_rlocs, axis = 1)  # copy endpoints
        # grad_n = np.append(grad_n_almost, np.array([grad_n_mid[:, -1]]).T, axis = 1)  # copy endpoints
        # aLne = a_minor * grad_n / dens
        
        
        # Plot density and density gradient
        
        fig = plt.figure(figsize = fig_size)
        if plot_eq:
            dens_sym_mid = dens_sym[:, :-1] + np.diff(dens_sym, axis = 1) / 2.
            aLne_eq_mid = a_minor * np.abs((np.diff(dens_sym, axis = 1) / np.diff(self.r_locs))) / dens_sym_mid
            # aLne_eq = np.append(aLne_eq_mid, np.array([aLne_eq_mid[:, -1]]).T, axis = 1)
            
            plt.pcolormesh(tor_grid[:,:-1], y_grid_orig_mid, aLne_eq_mid, cmap = cmap_d)
            
        else:
            plt.pcolormesh(tor_grid[:,:-1], y_grid_mid, aLne_mid, cmap = cmap_d)
                
        plt.colorbar()
        if clim_max is not None:
            plt.clim([0, clim_max])
        
        # Plot contours of ne (and Te if wanted)
        if plot_Te:
            if plot_eq:
                ct_Te = plt.contour(tor_grid, y_grid_orig, -Te_sym, nTe, colors = 'w')  # negative makes it dashed
            else:
                ct_Te = plt.contour(tor_grid, y_grid, -Te, nTe, colors = 'w') # negative makes it dashed
                
        if plot_psi_contours:
            if R_not_psi:
                if plot_eq:
                    # for p in psi_conts:
                    #     R_dash = f_psi_to_R_eq(p)
                    #     plt.plot([0,359], [R_dash, R_dash], '--w', lw = 1)
                    ct_psi = plt.contour(tor_grid, y_grid, -psi_grid_orig, colors = 'w', levels = -psi_conts)
                else:
                    # for p_ind, p in enumerate(psi_conts):
                    #     plt.plot(phi, R_dash_3D[:, p_ind], '--w', lw=1)
                    ct_psi = plt.contour(tor_grid, y_grid, -psi_grid, colors = 'w', levels = -psi_conts)

                # if dpsi_cont < 0.01:
                #     plt.clabel(ct_psi, fontsize = 12, fmt='%1.3f', inline = 0) #, color = 'w'
                # else:
                #     plt.clabel(ct_psi, fontsize = 12, fmt='%1.2f', inline = 0) #, color = 'w'
                            
            else:
                for p in psi_conts:
                    plt.plot([0,359], [p, p], '--w', lw = 1)  # No need to label, already on y-axis
            
        if plot_eq:
            ct_dens = plt.contour(tor_grid, y_grid_orig, dens_sym, ndens, colors = 'r')
        else:
            ct_dens = plt.contour(tor_grid, y_grid, dens, ndens, colors = 'r')

        if plot_diagnostics:
            plt.plot([145, 145], y_range, '--c', lw = 3)
            plt.plot([255, 255], y_range, '--c', lw = 3)
            if plot_labels:
                text_y_loc = y_range[0] + (y_range[1] - y_range[0]) / y_range[1] * 0.2
                plt.text(145, text_y_loc, 'BES', fontsize = 20, color = 'c', fontweight = 'bold',
                         horizontalalignment = 'right', verticalalignment = 'bottom', rotation='vertical')
                plt.text(255, text_y_loc, 'Reflectometer', fontsize = 20, color = 'c', fontweight = 'bold',
                         horizontalalignment = 'right', verticalalignment = 'bottom', rotation='vertical')
        
        if R_not_psi:
            # plt.yticks(np.arange(np.around(y_LCFS, decimals = 2), y_range[0] - dR/2, -dR))
            plt.yticks(np.arange(np.ceil(y_range[0]/dR)*dR, np.floor(y_range[1]/dR)*dR + dR / 2, dR))
            if plot_psi_contours or plot_Te:
                plt.plot(phi, R_LCFS_3D, '--w', lw=3)
        else:
            plt.plot([0, 359], [1, 1], '--w', lw=3)
            plt.yticks(np.arange(1., psi_range[0] - dpsi_label/2, -dpsi_label))

        plt.ylim(y_range)
        plt.xlim([0, 359])
        
        if plot_labels:
            plt.xlabel('$\\phi$ (deg)', fontsize = 20)
            if R_not_psi:
                plt.ylabel('R (m)', fontsize = 20)
            else:
                plt.ylabel('$\\psi_N$', fontsize = 20)
                
        plt.xticks(np.arange(0, 360, 60))
        plt.tight_layout()
        
        plt.tick_params(which = 'both', width = lw_border)
        plt.tick_params(which = 'major', length = 7)
        plt.tick_params(which = 'minor', length = 4)
        for axis in ['top', 'bottom', 'left', 'right']:
            fig.axes[0].spines[axis].set_linewidth(lw_border)
        
        # Now plot a/LTe for comparison
        
        if plot_aLTe:
            
            fig = plt.figure(figsize = fig_size)
            if plot_eq:
                Te_sym_mid = Te_sym[:, :-1] + np.diff(Te_sym, axis = 1) / 2.
                aLTe_eq_mid = a_minor * np.abs((np.diff(Te_sym, axis = 1) / np.diff(self.r_locs))) / Te_sym_mid
                # aLTe_eq = np.append(aLTe_eq_mid, np.array([aLTe_eq_mid[:, -1]]).T, axis = 1)
                plt.pcolormesh(tor_grid[:, :-1], y_grid_orig_mid, aLTe_eq_mid, cmap = cmap_t)
            else:
                Te_mid = Te[:, :-1] + np.diff(Te, axis = 1) / 2.
                aLTe_mid = a_minor * np.abs((np.diff(Te, axis = 1) / np.diff(self.r_locs))) / Te_mid
                # aLTe = np.append(aLTe_short, np.array([aLTe_short[:, -1]]).T, axis = 1)
                plt.pcolormesh(tor_grid[:, :-1], y_grid_mid, aLTe_mid, cmap = cmap_t)
            plt.colorbar()
            if clim_max is not None:
                plt.clim([0, clim_max])
                
            if plot_Te:
                if plot_eq:
                    ct_Te = plt.contour(tor_grid, y_grid_orig, -Te_sym, nTe, colors = 'w')  # negative makes it dashed
                else:
                    ct_Te = plt.contour(tor_grid, y_grid, -Te, nTe, colors = 'w')  # negative makes it dashed

            if R_not_psi:
                plt.yticks(np.arange(np.around(y_LCFS, decimals = 2), y_range[0] - dR / 2, -dR))
                plt.plot(phi, R_LCFS_3D, '--w', lw = 3)
            else:
                plt.plot([0, 359], [1, 1], '--w', lw = 3)
                plt.yticks(np.arange(1., psi_range[0] - dpsi_label / 2, -dpsi_label))
                
            plt.ylim(y_range)
            plt.xlim([0, 359])

            if plot_diagnostics:
                plt.plot([145, 145], y_range, '--w', lw = 2)
                plt.plot([255, 255], y_range, '--w', lw = 2)
                if plot_labels:
                    text_y_loc = y_range[0] + (y_range[1] - y_range[0]) / y_range[1] * 0.2
                    plt.text(145, text_y_loc, 'BES', fontsize = 20, color = 'w', fontweight = 'bold',
                             horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
                    plt.text(255, text_y_loc, 'Reflectometer', fontsize = 20, color = 'w', fontweight = 'bold',
                             horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
            if plot_labels:
                plt.xlabel('$\\phi$ (deg)', fontsize = 20)
                plt.ylabel('$\\psi_N$', fontsize = 20)
            plt.xticks(np.arange(0, 360, 60))
            plt.tight_layout()
            
        if plot_Ti:
            
            fig = plt.figure(figsize = fig_size)

            Ti_sym, _ = np.meshgrid(self.eqs[self.ti_ind, :], phi)  # toroidal locs X surfaces

            Ti = Ti_sym + perts_tot[:, self.ti_ind, :]  # total Ti (equilib + perturbed)
            Ti[np.where(Ti < Te_thresh)] = Te_thresh  # don't plot the negative Ti contours in the SOL

            if plot_eq:
                plt.pcolormesh(tor_grid, y_grid_orig, Ti_sym, cmap = cmap_t)
            else:
                plt.pcolormesh(tor_grid, y_grid, Ti, cmap = cmap_t)
                
            plt.colorbar()
            if clim_max is not None:
                plt.clim([0, clim_max])
                
            if plot_Te:
                if plot_eq:
                    ct_Te = plt.contour(tor_grid, y_grid_orig, -Te_sym, nTe, colors = 'w')  # negative makes it dashed
                else:
                    ct_Te = plt.contour(tor_grid, y_grid, -Te, nTe, colors = 'w')  # negative makes it dashed
                    
            # Overplot contours
            if plot_eq:
                ct_dens = plt.contour(tor_grid, y_grid_orig, Ti_sym, ndens, colors = 'r')
            else:
                ct_dens = plt.contour(tor_grid, y_grid, Ti, ndens, colors = 'r')

            if R_not_psi:
                plt.yticks(np.arange(np.around(y_LCFS, decimals = 2), y_range[0] - dR / 2, -dR))
                plt.plot(phi, R_LCFS_3D, '--w', lw = 3)
            else:
                plt.plot([0, 359], [1, 1], '--w', lw = 3)
                plt.yticks(np.arange(1., psi_range[0] - dpsi_label / 2, -dpsi_label))
                
            plt.ylim(y_range)
            plt.xlim([0, 359])

            if plot_diagnostics:
                plt.plot([145, 145], y_range, '--w', lw = 2)
                plt.plot([255, 255], y_range, '--w', lw = 2)
                if plot_labels:
                    text_y_loc = y_range[0] + (y_range[1] - y_range[0]) / y_range[1] * 0.2
                    plt.text(145, text_y_loc, 'BES', fontsize = 20, color = 'w', fontweight = 'bold',
                             horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
                    plt.text(255, text_y_loc, 'Reflectometer', fontsize = 20, color = 'w', fontweight = 'bold',
                             horizontalalignment = 'right', verticalalignment = 'bottom', rotation = 'vertical')
            if plot_labels:
                plt.xlabel('$\\phi$ (deg)', fontsize = 20)
                plt.ylabel('$\\psi_N$', fontsize = 20)
            plt.xticks(np.arange(0, 360, 60))
            plt.tight_layout()
            
        plt.show(block = False)

   
# -----------------------------------------------------------------------

class _M3DC1_dens_xsection(object):
    def __init__(self, sav_fileloc = '/home/wilcoxr/m3dc1/profile_savs/dens_xsection.sav',
                 linfac = -2.5 * 4 / np.pi):
        """
        I'm only saving this in case I have trouble making class M3DC1_xsection backwards compatible

        Requires previously saved IDL .sav file containing profile variables
        I generate this data with this (on venus): /u/wilcoxr/idl/load_m3dc1_single_n_fullxsection.pro

        ** This is only OK to do simply because I'm using the toroidal location at the peak of the density perturbation!
           Any other location requires more processing

        For original plot_dens_xsection: sav_fileloc = '/home/wilcoxr/m3dc1/profile_savs/dens_xsection.sav',
                                         linfac = -2.5*4/np.pi
        """
    
        idlsave_dict = si.readsav(sav_fileloc, python_dict = True, verbose = False)
    
        # fields = ['te','den','ti','pe','p'] (usually ordered like this in arrays)
    
    
        self.fields = idlsave_dict[
            'fields']  # fields that are in the IDL sav file, and which index corresponds to which
        self.r = idlsave_dict['r']  # grid locations (401)
        self.z = idlsave_dict['z']
        self.dens_pert = linfac * np.real(
            idlsave_dict['res'])  # Not correct unless just looking at the peak of the dens pert

    # -----------------------------------------------------------------------

    def plot_dens_xsection(self, cmap = 'bwr', no_labels = False, plot_wall = True,
                           plot_Icoils = True, clims = 3, norm = 1e18):
        """

        cmap
        no_labels
        plot_wall
        plot_Icoils
        clims        Should be given as a single number, and the color map will be symmetric about that
        """
    
        dens = self.dens_pert[:, :, 0] / norm
    
        f = plt.figure()
    
        plt.pcolormesh(self.r, self.z, dens, cmap = cmap)
    
        if plot_wall:
            R_wall, Z_wall = get_wall()
            plt.plot(R_wall, Z_wall, '-k', zorder = 1)
    
        if plot_Icoils:
            R_IL, Z_IL, R_IU, Z_IU = get_Icoils()
            plt.plot(R_IL, Z_IL, '-k', lw = 4, zorder = 2)
            plt.plot(R_IU, Z_IU, '-k', lw = 4, zorder = 2)
    
        if clims is None:
            clims = np.max(np.abs(dens))
    
        plt.clim([-clims, clims])
        cb = plt.colorbar()
        cb.set_ticks(range(-clims, clims + 1))
    
        plt.axis('equal')
        f.axes[0].set_xlim(min(self.r) - 0.1, max(self.r) + 0.1)
        f.axes[0].set_ylim(min(self.z) - 0.12, max(self.z) + 0.04)
        if no_labels:
            plt.xticks([1, 1.5, 2, 2.5], ['', '', '', ''])
            plt.yticks([-1, -0.5, 0, 0.5, 1], ['', '', '', '', ''])
        else:
            plt.xticks([1, 1.5, 2, 2.5])
            plt.yticks([-1, -0.5, 0, 0.5, 1])
            plt.xlabel('R (m)')
            plt.ylabel('Z (m)')
    
        plt.tight_layout()

    # -----------------------------------------------------------------------

    def plot_Te_perts(self, tor = 45, clims = None, no_labels = False):
        f = plt.figure()
        cb = plt.colorbar()
        if clims is not None:
            plt.clim([-clims, clims])
            cb.set_ticks(range(-clims, clims + 1))
    
        plt.axis('equal')
        f.axes[0].set_xlim(min(self.r) - 0.1, max(self.r) + 0.1)
        f.axes[0].set_ylim(min(self.z) - 0.12, max(self.z) + 0.04)
        if no_labels:
            plt.xticks([1, 1.5, 2, 2.5], ['', '', '', ''])
            plt.yticks([-1, -0.5, 0, 0.5, 1], ['', '', '', '', ''])
        else:
            plt.xticks([1, 1.5, 2, 2.5])
            plt.yticks([-1, -0.5, 0, 0.5, 1])
            plt.xlabel('R (m)')
            plt.ylabel('Z (m)')
    
        plt.tight_layout()
        
# # -----------------------------------------------------------------------
#
# def plot_dens_turb(max_amp = 0.005, n = 3, psi_range = [0.95, 1.01], dpsi = 0.01, ped_top = 0.97, surf_amp = 0.0005,
#                    fig_size = (8, 7.5)):
#     npts_rad = 200
#     npts_tor = 500
#     phi = np.linspace(0, 360, npts_tor)
#     phi_rad = np.deg2rad(phi)
#
#     surf_locs = np.arange(1., psi_range[0], -dpsi)
#
#     dens = np.ones([len(surf_locs)])
#
#     # f, (ax1, ax2) = plt.subplots(2, 1, figsize = fig_size)
#     fig = plt.figure(figsize = fig_size)
#
#     for s in surf_locs:
#         surf_contour = s + surf_amp * np.sin(n * phi_rad)
#
#         amp = np.max([0, max_amp * (s - ped_top) / (1 - ped_top)])
#         dens_contour = s + surf_amp * np.sin(n * phi_rad) + amp * np.sin(n * phi_rad)
#         turb_contour = s + surf_amp * np.sin(n * phi_rad) - amp * np.sin(n * phi_rad)
#
#         if s == 1:
#             lw = 2
#         else:
#             lw = 1
#         plt.plot(surf_contour, phi, '--k', lw = lw)
#
#         plt.plot(dens_contour, phi, 'r')
#         plt.plot(turb_contour, phi, '-.g')
#
#     plt.xlim(psi_range)
#     plt.ylim([0, 360])
#     plt.ylabel('$\\phi$ (deg)', fontsize = 14)
#     plt.yticks(np.arange(0, 361, 60))
#     plt.xticks(surf_locs)
#
#     # plt.tight_layout()
    
#-----------------------------------------------------------------------

def __plot_dens_turb_contour(max_dens_amp = 3, n = 3, psi_range = [0.924, 1.009], dpsi=0.01, ped_top = 0.97,
                           sep_dens = 10, offset_deg = 0,
                           tanh_wid = 0.03, surf_amp = 0.0015, ndens = 20, max_turb_amp = 0.009, fig_size=(13, 8),
                           plot_diagnostics = True, plot_profiles = False, plot_turb_contours = False, plot_labels = True):
    
    """
    Plots cartoon of how this might look (now should use plot_surf_contours_OMP instead)
    
    Inputs:
    
    max_dens_amp         Maximum amplitude of the perturbed density (at separatrix, decreases linearly to the top of
                         the pedestal from there)
    n                    Toroidal mode number of applied 3D field
    psi_range
    dpsi                 Differential flux surface spacing for plotting
    ped_top              Location of pedestal top in normalized flux
    sep_dens             Separatrix density, in 10^18 m^-3
    offset_deg           Should be ~5 degrees, but then y-axis doesn't line up, so just put 0
    tanh_wid             Width of the tanh function for the density
    surf_amp             Amplitude of the surface deformation, in normalized flux
    ndens                Number of density contours
    max_turb_amp
    fig_size
    plot_diags
    plot_profiles
    plot_turb_contours

    """

    npts_rad = 200
    npts_tor = 500
    phi = np.linspace(0, 360, npts_tor)
    phi_rad = np.deg2rad(phi)
    rad_locs = np.linspace(psi_range[0], psi_range[1], npts_rad)
    
    surf_locs = np.arange(1., psi_range[0], -dpsi)
    
    dens_lin = 6 / 0.15 # (e18 m^-3 / psi_n)

    dens_profile = sep_dens + np.maximum(0, 1 - rad_locs) * dens_lin + sep_dens * np.tanh((1 - rad_locs) / tanh_wid)
    
    if plot_profiles:
        dens_profile_max = dens_profile + np.minimum(max_dens_amp,
                                                     np.maximum(0, max_dens_amp * (rad_locs - ped_top) / (1 - ped_top)))
        dens_profile_min = dens_profile - np.minimum(max_dens_amp,
                                                     np.maximum(0, max_dens_amp * (rad_locs - ped_top) / (1 - ped_top)))
    
        plt.figure()
        plt.plot([1, 1], [0, 25], '--k')
        plt.plot(rad_locs, dens_profile, '--k', lw = 2)
        plt.plot(rad_locs, dens_profile_max, 'r', lw = 2)
        plt.plot(rad_locs, dens_profile_min, 'b', lw = 2)
        plt.xlim(psi_range)
        plt.ylabel('n$_e$ (10$^{18}$ m$^{-3}$)')
        plt.xlabel('$\\psi_N$')
    
    tor_grid, rad_grid = np.meshgrid(phi_rad, rad_locs)  # surfaces X toroidal locs
    tor_grid_deg = tor_grid * 180 / np.pi
    dens = sep_dens + np.maximum(0, 1 - rad_grid) * dens_lin + sep_dens * np.tanh((1 - rad_grid) / tanh_wid)

    dens_amp = np.minimum(max_dens_amp, np.maximum(0, max_dens_amp * (rad_grid - ped_top) / (1 - ped_top)))
    dens_surf_amp = -surf_amp * np.diff(dens_profile) / np.diff(rad_locs)
    _, dens_surf_amp_grid = np.meshgrid(phi_rad, np.append(dens_surf_amp, dens_surf_amp[-1]))
    
    dens += dens_surf_amp_grid * np.sin(n * (tor_grid + offset_deg * np.pi / 180))
    dens += dens_amp * np.sin(n * (tor_grid + offset_deg * np.pi / 180))
    

    fig = plt.figure(figsize = fig_size)
    ct = plt.contour(tor_grid_deg, rad_grid, dens, ndens, colors = 'r')
    
    ds = 0.001  # to separate turbulence contours and magnetic surfaces
    for s in surf_locs:
        surf_contour = s + surf_amp * np.sin(n * (phi_rad + offset_deg * np.pi / 180))
        
        if plot_turb_contours:
            turb_amp = np.max([0, max_turb_amp * (s - ped_top) / (1 - ped_top)])
            turb_contour = s + ds + surf_amp * np.sin(n * (phi_rad + offset_deg * np.pi / 180)) \
                           - turb_amp * np.sin(n * (phi_rad + offset_deg * np.pi / 180))
            
            plt.plot(phi, turb_contour , ':b', lw = 2)
        
        if s == 1:
            lw = 2
        else:
            lw=1
        plt.plot(phi, surf_contour, '--k', lw = lw)
        
    if plot_diagnostics:
        plt.text(145, 1 + surf_amp, '$\\downarrow$', fontsize = 30, color='k', fontweight = 'bold',
                 horizontalalignment='center', verticalalignment='bottom')
        plt.text(255, 1 + surf_amp, '$\\downarrow$', fontsize = 30, color='k', fontweight = 'bold',
                 horizontalalignment='center', verticalalignment='bottom')
        if plot_labels:
            plt.text(152, 1 + surf_amp, 'BES', fontsize = 20, color='k', fontweight = 'bold',
                     horizontalalignment='left', verticalalignment='bottom')
            plt.text(262, 1 + surf_amp, 'Reflectometer', fontsize = 20, color='k', fontweight = 'bold',
                     horizontalalignment='left', verticalalignment='bottom')

        
    plt.ylim(psi_range)
    plt.xlim([0, 360])
    if plot_labels:
        plt.xlabel('$\\phi$ (deg)', fontsize = 20)
        plt.ylabel('$\\psi_N$', fontsize = 20)
    plt.xticks(np.arange(0, 361, 60))
    plt.yticks(surf_locs)
    plt.tight_layout()

# -----------------------------------------------------------------------
    
def make_segmented_cmap():
    """
    From Stack Overflow:
    http://stackoverflow.com/questions/23712207/cyclic-colormap-without-visual-distortions-for-use-in-phase-angle-plots
    """
    white = '#ffffff'
    black = '#000000'
    red = '#ff0000'
    blue = '#0000ff'
    yellow = '#ffff00'
    anglemap = col.LinearSegmentedColormap.from_list(
        'anglemap', [black, red, yellow, blue, black], N=256, gamma=1)
    return anglemap

# --- get_wall and get_Icoils ------------------------------------------

def get_wall():
    Rwall = [1.41903, 1.28008, 1.27999, 1.24803, 1.22784, 1.20913, 1.19011, 1.16185,
             1.11593, 1.04209, 1.02923, 1.00088, 1.01777, 1.01779, 1.01619, 1.01627,
             1.01608, 1.15285, 1.19734, 1.24185, 1.28635, 1.33085, 1.37535, 1.41985,
             1.41985, 1.372, 1.372, 1.372, 1.57, 1.76801, 1.76801, 1.78603,
             2.13401, 2.37704, 2.35342, 2.35116, 2.354, 2.35082, 2.3523, 2.37704,
             2.1282, 2.0699, 1.78499, 1.647, 1.60799, 1.37198, 1.3721, 1.41897, 1.41903]
    Zwall = [1.34783, 1.34773, 1.33104, 1.2539, 1.22143, 1.20165, 1.18772, 1.1759,
             1.16606, 1.16223, 1.21657, 1.21721, 1.13839, 1.00132, 1.00132, -1.21079,
             -1.22884, -1.3664, -1.3664, -1.3664, -1.3664, -1.3664, -1.3664, -1.3664,
             -1.329, -1.329, -1.25, -1.25, -1.25, -1.25, -1.21104, -1.17428,
             -0.973024, -0.389023, -0.400247, -0.337078, 0., 0.204946, 0.400178, 0.389023,
             0.993424, 1.03975, 1.07688, 1.07675, 1.09525, 1.29208, 1.30954, 1.31017, 1.34783]
    return Rwall, Zwall

def get_Icoils():
    R_IL = [2.16408, 2.37363]
    Z_IL = [-1.0119614, -0.5033518]
    
    R_IU = [2.37363, 2.164156]
    Z_IU = [0.5033518, 1.01196]
    
    return R_IL, Z_IL, R_IU, Z_IU

def get_pellet_traj(port = 'HFS_upper'):
    if port.upper() == 'HFS_UPPER':
        R_traj = [1.016, 1.6]
        Z_traj = [0.721, 0]
    elif port.upper() == 'HFS_MID':
        R_traj = [1.016, 1.6]
        Z_traj = [0.280, 0]
    elif  port.upper() == 'OUTBOARD':
        R_traj = [2.365, 1.15]
        Z_traj = [0.089, 0.089]
    elif  port.upper() == 'SPI':
        R_traj = [2.284, 1.016]
        Z_traj = [0.684, -0.690]
    elif  port.upper() == 'R2':
        R_traj = [2.246, 1.016]
        Z_traj = [-1.474, -0.007]
    elif  port.upper() == 'V1':
        R_traj = [1.479, 1.479]
        Z_traj = [1.203, 0.0]
    elif  port.upper() == 'V3':
        R_traj = [2.098, 2.098]
        Z_traj = [1.017, 0.0]
    else:
        print 'Requested port not available: ' + port
        return
        
    return R_traj, Z_traj

