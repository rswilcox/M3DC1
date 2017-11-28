import numpy as np
import os
import matplotlib.pyplot as plt


# plt.rcParams.update({'font.family': 'Arial'})
plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'figure.facecolor': 'w'})
plt.rcParams.update({'mathtext.default': 'regular'})


OMFIT_M3DC1_folder = '/fusion/projects/results/m3dc1_results/wilcoxr/data/OMFIT/DIII-D/'


#-----------------------------------------------------------------------

class Profile(object):
    def __init__(self, fileloc, verbose = False):
        """
        Object of an M3D-C1 profile file, which can be plotted or modified
        """
        
        # with open(fileloc, 'r') as f:
        #     lines = f.readlines()
        
        
        profile = np.loadtxt(fileloc)  # shape = (n, 2)

        self.psin = profile[:,0]
        self.values = profile[:,1]
        
        self.profile_name = get_profile_name(fileloc)
    
    # -------------------------------------------------------------------
    
    def plot_profile(self):
        plt.figure()
        plt.plot(self.psin, self.values, 'k', lw=2)
        plt.xlabel('$\\psi_n$')
        plt.ylabel(self.profile_name)
        plt.tight_layout()

    # ------------------------------------------------------------------
    
    def add_gauss_edge(self, gauss_peak, psin_edge = 1.1, gauss_wid = 0.001,
                       plot_change = True, saveit = False):

        old_psin = self.psin.copy()
        old_values = self.values.copy()
        
        self.extend_range(psin_edge, plotit = False, saveit = saveit)
            
        gaussian_add = gauss(self.psin, psin_edge, gheight = gauss_peak, gwidth = gauss_wid)
            
        if plot_change:
            plt.figure()
            plt.plot(old_psin,old_values, 'b')
            plt.plot(self.psin, self.values + gaussian_add, 'r')
            plt.xlim([0.9, psin_edge + 0.01])
            plt.xlabel('$\\psi_n$')
            plt.ylabel(self.profile_name)
            plt.tight_layout()
            
        if saveit: self.values += gaussian_add

    # ------------------------------------------------------------------
            
    def extend_linearly(self, psin_end = 1.02, saveit = False, plot_change = True):
        old_psin = self.psin.copy()
        old_values = self.values.copy()
    
        new_psin, new_values = self.extend_range(psin_end, plotit = False, saveit = False)
        
        end_gradient = (old_values[-1] - old_values[-2]) / (old_psin[-1] - old_psin[-2])
        
        outside_inds = (new_psin > np.max(old_psin))
        new_values[outside_inds] = old_values[-1] + \
                                   (new_psin[outside_inds] - old_psin[-1]) * end_gradient
                                   
        if plot_change:
            plt.figure()
            plt.plot(old_psin, old_values, 'b')
            plt.plot(old_psin[-1], old_values[-1], 'ob')
            plt.plot(new_psin, new_values, 'r')
            plt.plot(new_psin[-1], new_values[-1], 'or')
            plt.xlim([0.9, psin_end + 0.01])
            plt.xlabel('$\\psi_n$')
            plt.ylabel(self.profile_name)
            plt.tight_layout()
            
        if saveit:
            self.psin = new_psin
            self.values = new_values
            
            
    # ------------------------------------------------------------------
            
    def extend_range(self, psin_end, psin_from_which_to_extrap_psi = 0.97,
                     plotit = False, saveit=False, verbose = True):
        """
        Just extends the psin and values arrays out to a further value
        values are copied from the last one out
        
          psin_end    This value of psin will not be exceeded
        """
        dpsi = np.mean(np.diff(self.psin[self.psin > psin_from_which_to_extrap_psi]))

        if np.max(self.psin) + dpsi <= psin_end:
            new_psin = np.append(self.psin, np.arange(self.psin[-1] + dpsi, psin_end, dpsi))
            new_npoints = len(new_psin) - len(self.values)
            new_values = np.append(self.values, np.ones(new_npoints) * self.values[-1])
            
        else:
            if verbose: print "Range already extended that far"
            return self.psin.copy(), self.values.copy()
        
        if plotit:
            plt.figure()
            plt.plot(self.psin, self.values, 'b', lw=2)
            plt.plot(self.psin[-1], self.values[-1], 'ob', lw = 2)
            plt.plot(new_psin, new_values, 'r')
            plt.plot(new_psin[-1], new_values[-1], 'or')
            plt.xlim([0.9, psin_end + dpsi])
            plt.xlabel('$\\psi_n$')
            plt.ylabel(self.profile_name)
            plt.title('Extended range')
            plt.tight_layout()
        
        if saveit:
            self.psin = new_psin
            self.values = new_values
            
        return new_psin, new_values

    # ------------------------------------------------------------------
            
    def write_new_profile(self, new_fileloc):
        """
        This at least can reproduce data that it reads...
        """
        
        with open(new_fileloc, 'w') as f:
            for i in range(len(self.psin)):
                f.write('  {0:>8}   {1:>17}  \n'.format(self.psin[i], self.values[i]))
                
        
        
        
#-----------------------------------------------------------------------

def compare_profiles(fileloc_list, plot_final_point = True):
    line_colors = 'brkgcmy'
    plt.figure()
    
    for i, f in enumerate(fileloc_list):
        p = Profile(f)
        
        # path, filename = os.path.split(f)
        # efit_ind = path.index('efit')
        # run_ident = path[]
        
        plt.plot(p.psin, p.values, line_colors[i], lw=2, label = str(i+1))
        if plot_final_point:
            plt.plot(p.psin[-1], p.values[-1], 'o'+line_colors[i])
            
        
    # plt.xlim([0.9, psin_edge + dpsi])
    plt.xlabel('$\\psi_n$')
    plt.ylabel(p.profile_name)
    plt.legend()
    plt.tight_layout()
    

#-----------------------------------------------------------------------
        
def get_profile_name(fileloc):
    path, filename = os.path.split(fileloc)
    
    if filename[:8] != 'profile_':
        print "Warning, filename does not match convention, should start with 'profile_'"
        
    return filename[8:]

#-----------------------------------------------------------------------
    
def gauss(x, gloc, gheight, gwidth):
    return gheight * np.exp( -(x - gloc )**2 / (2*gwidth) )