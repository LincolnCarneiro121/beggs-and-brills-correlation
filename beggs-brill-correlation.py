from math import pi,log,sin,radians,log10, exp
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
    reservoir_pressure = read_float_input('Enter the reservoir pressure (psia): ')
    reservoir_temperature = read_float_input('Enter reservoir temperature (ºF): ')
    GOR = read_float_input('Enter GOR (scf/stb): ')
    Bo = read_float_input('Enter Bo (bbl/stb): ')
    oil_viscosity = read_float_input('Enter oil viscosity (cp): ')
    surface_tension = read_float_input('Enter surface tension (dina/cm²): ')
    gas_density = read_float_input('Enter gas density: ')
    oil_density = read_float_input('Enter oil density: ')
    gas_ratio = read_float_input('Enter gas ratio (ft³/s): ')
    liquid_ratio = read_float_input('Enter liquid ratio (ft³/s): ')
    tubing_area = read_float_input('Enter tubing area (ft²): ')
    theta = read_float_input('Enter the inclination angle (degrees): ')
    inclination = read_inclination_string_input('Descendent[d] or ascendent[a]? ')

    #reservoir_pressure = 1719.7
    #reservoir_temperature = 90
    #GOR = 947.5
    #Bo = 1.495
    #oil_viscosity = 0.5
    #surface_tension = 28
    #gas_density = 8.823
    #oil_density = 38.32
    #gas_ratio = 0.08855
    #liquid_ratio = 0.0466
    #tubing_area = 0.0217
    #theta = 90
    #inclination = 'a'

    # gas gravity
    gas_gravity = round(gas_density/16.02,4)
    print('gas_gravity: ', gas_gravity)

    # oil gravity
    oil_gravity = round(oil_density/62.4,4)
    print('oil_gravity: ', oil_gravity)

    # gravity (ft²/s)
    g = 32.174

    # gas viscosity (cp)
    Mg = 28.967*gas_gravity
    K1 = ((0.00094 + 2*10**(-6)*Mg)*(reservoir_temperature + 459.67)**1.5)/(209 + 19*Mg + 459.67 + reservoir_temperature)
    X = 3.5 + 986/(459.67 + reservoir_temperature) + 0.01*Mg
    Y = 2.4 - 0.2*X
    gas_viscosity = round(K1*exp(X*0.016018*gas_density**Y),6)

    print(Mg, ' / ', K1, ' / ', X, ' / ', Y)
    print('gas_viscosity: ', gas_viscosity, ' cp')

    # diameter
    ID = round(((4*tubing_area)/(pi))**(1/2),4)
    print('Internal Tubing Diameter: ', ID, ' ft')

    # surface velocity
    vsl = round(liquid_ratio/tubing_area,4)
    vsg = round(gas_ratio/tubing_area,4)
    vm = round(vsl + vsg,2)

    print('vsl, vsg, vm: ', vsl, ' / ', vsg, ' / ', vm, ' ft/s')

    # no slip hold up
    Cl = round(liquid_ratio/(gas_ratio + liquid_ratio),4)
    print('no slip: ', Cl)

    # Froude Number
    Froude = round((vm**2)/(g*ID),4)
    print('Froude: ', Froude)

    # mixture viscosity
    viscosity = round((oil_viscosity*Cl + (1-Cl)*gas_viscosity)*(6.72*10**(-4)),7)
    print('viscosity: ', viscosity, ' cp')


    Gt = round(vsl*oil_density + vsg*gas_density,4)

    # Reynolds number
    Re = round((ID*Gt)/viscosity,4)
    print('Reynolds: ', Re)

    # velocity number
    Nvl = round(1.938*vsl*(oil_density/(surface_tension))**(1/4),4)
    print('velocity number Nvl: ', Nvl)

    # flow pattern
    L1 = round(316*Cl**0.302,4)
    L2 = round(0.0009252*Cl**(-2.4684),4)
    L3 = round(0.1*Cl**(-1.4516),4)
    L4 = round(0.5*Cl**(-6.738),4)

    print('L1, L2, L3, L4: ', L1, '/', L2, '/', L3, '/', L4)

    global beta

    if (Cl < 0.01 and Froude < L1) or (Cl >= 0.01 and Froude < L2):
        print('SEGREGATED FLOW')
        EL_0 = round((0.98*Cl**0.4846)/(Froude**0.0173),4)
        if inclination == 'a':
            beta = round((1 - Cl)*log((0.011*Nvl**3.539)/((Cl**3.768)*Froude**1.614)),4)

    if (0.01 <= Cl < 0.4 and L3 < Froude <= L1) or (Cl >= 0.4 and L3 < Froude <= L4):
        print('INTERMITENT FLOW')
        EL_0 = round((0.845*Cl**0.5351)/(Froude**0.0173),4)
        if inclination == 'a':
            beta = round((1 - Cl)*log(((2.96*Cl**0.305)*Froude**0.0978)/(Nvl**0.4473)),4)


    if (Cl < 0.4 and Froude >= L4) or (Cl >= 0.4 and Froude > L4):
        print('DISTRIBUTED FLOW')
        EL_0 = round((1.065*Cl**0.5824)/(Froude**0.0609),4)
        if inclination == 'a':
            beta = 0

    if L2 < Froude < L3:
        print('TRANSITION FLOW')
        EL_0 = ((L3-Froude)/(L3-L2))*((0.98*Cl**0.4846)/(Froude**0.0173)) + (1-(L3-Froude)/(L3-L2))*((0.845*Cl**0.5351)/(Froude**0.0173))
        return
    
    if inclination == 'd':
        beta = round((1 - Cl)*log((4.7*Nvl**0.1244)/((Cl**0.3692)*Froude**0.5056)),4)

    # inclination correction factor
    B_theta = round(1 + beta*(sin(radians(1.8*theta)) - (1/3)*sin(radians(1.9*theta)**3)),4)

    print('Holdup (0): ', EL_0)
    print('beta or C: ', beta)
    print('inclination correction factor: ', B_theta)

    # liquid hold up
    EL_theta = round(B_theta*EL_0,4)

    print('liquid holdup: ', EL_theta)

    # mixture density
    mixture_density = round(oil_density*EL_theta + gas_density*(1-EL_theta),4)

    print('mixture_density: ', mixture_density, ' lb/ft³')

    # friction factor ratio
    y = round(Cl/(EL_theta**2),4)
    
    ln_y = round(log(y),4) 

    print('y: ', y)
    print('ln(y): ', ln_y)

    ft_fns = round(exp(ln_y/(-0.0523 + 3.182*ln_y - 0.8725*ln_y**2 +0.1853*ln_y**4)),4)

    print('ft/fns: ', ft_fns)

    # no slip friction factor
    fns = 1/(2*log10(Re/(4.5223*log10(Re) - 3.8215)))**2

    print('fns: ', fns)
    # two phase friction factor

    ft = fns*ft_fns

    print('ft: ', ft)

    # pressure gradient
    dl_dp = 144*(1 - (mixture_density*vm*vsg)/(g*reservoir_pressure))/(mixture_density*sin(radians(theta)) + (ft*Gt*vm)/(2*g*ID))
    print('pressure gradient: ', round(1/dl_dp,2), 'psi/ft')


if __name__ == '__main__':
    main()