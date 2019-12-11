

def ziegler_files():
    path= '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/'

    SS_n10 = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_-10_fluxes.csv'
    SS_n5  = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_-5_fluxes.csv'
    SS_0   = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_0_fluxes.csv'
    SS_p5   = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_+5_fluxes.csv'
    SS_p10  = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_+10_fluxes.csv'

    Ni_n10  = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_Ni_-10_fluxes.csv'
    Ni_n5 =  '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_Ni_-5_fluxes.csv'
    Ni_0 ='/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_Ni_0_fluxes.csv'
    Ni_p5 ='/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_Ni_+5_fluxes.csv'
    Ni_p10 = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_Ni_+10_fluxes.csv'

    Cu_n10 = path + 'E_foils_Cu_-10_fluxes.csv'
    Cu_n5 = path + 'E_foils_Cu_-5_fluxes.csv'
    Cu_0 = path + 'E_foils_Cu_0_fluxes.csv'
    Cu_p5 = path + 'E_foils_Cu_+5_fluxes.csv'
    Cu_p10 = path + 'E_foils_Cu_+10_fluxes.csv'

    Fe_n10 = path + 'E_foils_Cu_-10_fluxes.csv'
    Fe_n5 = path + 'E_foils_Cu_-5_fluxes.csv'
    Fe_0 = path + 'E_foils_Cu_0_fluxes.csv'
    Fe_p5 = path + 'E_foils_Cu_+5_fluxes.csv'
    Fe_p10 = path + 'E_foils_Cu_+10_fluxes.csv'


    files = [SS_n10, SS_n5, SS_0, SS_p5, SS_p10, Ni_n10, Ni_n5, Ni_0, Ni_p5, Ni_p10, Cu_n10, Cu_n5, Cu_0, Cu_p5, Cu_p10, Fe_n10, Fe_n5, Fe_0, Fe_p5, Fe_p10,]
    names = ['-SS-10%', '-SS-5%', '-SS0%', '-SS+5%', '-SS+10%', '-Ni-10%', '-Ni-5%', '-Ni0%', '-Ni+5%', '-Ni+10%', '-Cu-10%', '-Cu-5%', '-Cu0%', '-Cu+5%', '-Cu+10%', '-Fe-10%', '-Fe-5%', '-Fe0%', '-Fe+5%', '-Fe+10%']
    return files, names
