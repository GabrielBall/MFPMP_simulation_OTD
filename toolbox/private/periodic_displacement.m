function displacement = periodic_displacement(displacement,domain)

domainWidth = domain.upper-domain.lower;
displacement = mod(displacement+0.5.*domainWidth,domainWidth)-0.5.*domainWidth;

end