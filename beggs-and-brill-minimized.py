from PIL import Image
from math import pi,log,sin,radians,log10
# Beggs and Brills Correlation

# General field info

def read_float_input(prompt):
    while True:
        try:
            value = float(input(prompt))
            return value
        except ValueError:
            print("Erro de digitação. Por favor, tente novamente.")

def read_inclination_string_input(prompt):
    while True:
        try:
            value = input(prompt)
            
            if value not in ['a','d']:
                raise ValueError
            
            return value
        except ValueError:
            print("Digite a ou d")

def main():
    oil_rate = read_float_input('Enter oil rate (stb/d): ')
    GOR = read_float_input('Enter Gas Oil ratio [GOR] (scf/stb): ')
    gas_gravity = read_float_input('Enter gas gravity (relative): ')
    API = read_float_input('Enter oil ºAPI: ')
    ID = read_float_input('Enter internal tubing diameter (in): ')
    BHP = read_float_input('Enter bottom hole pressure (psi): ')
    BHT = read_float_input('Enter bottom hole temperature (ºF): ')
    WHT = read_float_input('Enter well head temperature (ºF): ')
    theta = read_float_input('Enter the inclination angle (degrees): ')
    inclination = read_inclination_string_input('Descendent[d] or ascendent[a]? ')

    # State propierties

    oil_gravity = 141.5/(API + 131.5) 
    liquid_density = 62.4*oil_gravity # lb/ft³
    gas_density = 16.02*gas_gravity # lb/ft³

    print('oil gravity: ', oil_gravity)
    print('liquid density: ', liquid_density, ' lb/ft³')
    print('gas density: ', gas_density, ' lb/ft³')

    # viscosity (Kartoatmodjo & Schmidt)

    viscosity = (16*10**8)*(BHT**(-2.8177)*log10(API)**(5.7526*log10(BHT)-26.9718))

    print('oil viscosity: ', viscosity, ' cp')

    # surface tension of oil (Beggs)

    tension_68 = 39 - 0.2571*API
    tension_100 = 37.5 - 0.2571*API
    intermediate_surface_tension = tension_68 - ((WHT - 68)*(tension_68 - tension_100))/32

    if WHT <= 68:
        surface_tension = tension_68
    elif 68 <= WHT <= 100:
        surface_tension = intermediate_surface_tension
    else:
        surface_tension = tension_100

    print('surface tension: ', surface_tension, ' dina/cm²')


    Ppc = 706 + 51.7*gas_gravity - 11.1*gas_gravity**2 # psi
    Tpc = 187 + 330*gas_gravity - 71.5*gas_gravity**2 # R
    Ppr = BHP/Ppc
    Tpr = (BHT+459.67)/Tpc
    
    print('PPR: ', Ppr)
    print('TPR: ', Tpr)

    im = Image.open("fator-compressibilidade.PNG") 
    im.show()

    compressibility_factor = read_float_input('Enter compressibility factor: ')

    Bg = 0.0283*compressibility_factor*(BHT/BHP)

    print('Bg: ', Bg)
    
    # no-slip hold-up


    GOR_conv = GOR*0.1781 # scf/scf
    ID_conv = ID/12 # ft
    oil_rate_conv = 0.000065*oil_rate # ft³/s
    gas_rate_conv = GOR_conv*oil_rate_conv # ft³/s

    vsl = (4*oil_rate_conv)/(pi*ID_conv**2)
    vsg = (4*gas_rate_conv*Bg)/(pi*ID_conv**2)

    print('vsl: ',vsl, ' ft/s')
    print('vsg: ', vsg, ' ft/s')

    vm = vsl + vsg
    Cl = vsl/vm

    print('no-slip: ', Cl)

    L1 = 316*Cl**0.302
    L2 = 0.0009252*Cl**(-2.4684)
    L3 = 0.1*Cl**(-1.4516)
    L4 = 0.5*Cl**(-6.738)

    Froude = (vm**2)/(32.2*ID_conv)

    print('L1: ',L1)
    print('L2: ',L2)
    print('L3: ',L3)
    print('L4: ',L4)
    print('Froude: ', Froude)

    # liquid velocity number

    Nvl = 1.938*vsl*(liquid_density/(32.2*surface_tension))**(1/4)

    # Liquid hold-ups

    if (Cl < 0.01 and Froude < L1) or (Cl >= 0.01 and Froude < L2):
        print('SEGREGATED FLOW')
        EL_0 = (0.98*Cl**0.4846)/(Froude**0.0173)
        if inclination == 'a':
            beta = (1 - Cl)*log((0.011*Nvl**3.539)/((Cl**3.768)*Froude**1.614))

    if (0.01 <= Cl < 0.4 and L3 < Froude <= L1) or (Cl >= 0.4 and L3 < Froude <= L4):
        print('INTERMITENT FLOW')
        EL_0 = (0.845*Cl**0.5351)/(Froude**0.0173)
        if inclination == 'a':
            beta = (1 - Cl)*log(((2.96*Cl**0.305)*Froude**0.0978)/(Nvl**0.4473))


    if (Cl < 0.4 and Froude >= L4) or (Cl >= 0.4 and Froude > L4):
        print('DISTRIBUTED FLOW')
        EL_0 = (1.065*Cl**0.5824)/(Froude**0.0609)
        if inclination == 'a':
            beta = 0
    
    if inclination == 'd':
        beta = (1 - Cl)*log((4.7*Nvl**0.1244)/((Cl**0.3692)*Froude**0.5056))

    if L2 < Froude < L3:
        print('TRANSITION FLOW')
        EL_0 = ((L3-Froude)/(L3-L2))*((0.98*Cl**0.4846)/(Froude**0.0173)) + (1-(L3-Froude)/(L3-L2))*((0.845*Cl**0.5351)/(Froude**0.0173))

    if beta < 0:
        beta = 0

    B_theta = 1 + beta*(sin(radians(1.8*theta)) - (1/3)*sin(radians(1.9*theta)**3))
    
    EL_theta = B_theta*EL_0

    mixture_density = liquid_density*EL_theta + gas_density*(1-EL_theta)

    # hydrostatic loss

    hydrostatic_loss = mixture_density*sin(radians(theta))/144

    print('Corrected Hold-up: ', EL_theta)
    print('mixture density: ', mixture_density, ' lb/ft³')
    print('hydrostatic loss: ', hydrostatic_loss)

    # friction factor

    Re = 1488*(liquid_density*vsl*ID_conv)/viscosity
    stainless_steel_rugosity = 0.00000656
    e_ID = stainless_steel_rugosity/ID_conv

    print('Reynolds: ', Re)

    A = (2.457*log(1/((7/Re)**0.9 + 0.27*e_ID)))**16
    B = (37530/Re)**16

    f = 8*((8/Re)**12 + 1/((A + B)**(3/2)))**(1/12)

    no_slip_density = Cl*liquid_density + (1 - Cl)*gas_density

    friction_loss = (2*f*no_slip_density*vm**2)/(144*32.2*ID_conv)

    print('friction loss: ', friction_loss)

    print('total loss: ',  friction_loss + hydrostatic_loss, ' psi/ft')

if __name__ == '__main__':
    main()