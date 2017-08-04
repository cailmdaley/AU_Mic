from astrocail import fitting 
import emcee_fixed_starflux as fs

def rename(Chain):
    Chain.chain['starflux'] *= 1e6
    Chain.chain['d_r'] += Chain.chain['r_in']
    Chain.chain.rename(columns={
        'm_disk'       : '$\log$ m_disk',
        'r_in'         : 'r$_{in}$',
        'd_r'          : 'r$_{out}$',
        'starflux' : 'Starflux ($\mu$Jy)'}, inplace=True)
        
    Chain.clean_chain['starflux'] *= 1e6
    Chain.clean_chain['d_r'] += Chain.clean_chain['r_in']
    Chain.chain.rename(columns={
        'm_disk'       : '$\log$ m_disk',
        'r_in'         : 'r$_{in}$',
        'd_r'          : 'r$_{out}$',
        'starflux' : 'Starflux ($\mu$Jy)'}, inplace=True)
        
        
run7 = fitting.MCMCrun('run7', nwalkers=18)
reload(fitting)
reload(fs)
fs.make_best_fit(run7)
        





        
