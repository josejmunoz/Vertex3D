function [Energies,Energy]=UpdateEnergy(Set,incr,Energy)
% Update different energy values
if incr==0 
    Energies.S=zeros(Set.Nincr,1); % Surface
    Energies.V=zeros(Set.Nincr,1); % Volume
    Energies.F=zeros(Set.Nincr,1); % Friction/viscous
    Energies.b=zeros(Set.Nincr,1); % bending
    Energies.B=zeros(Set.Nincr,1); % Barrier
    Energies.C=zeros(Set.Nincr,1); % Contractility
    Energies.I=zeros(Set.Nincr,1); % Contractility at junctions
    Energies.Sub=zeros(Set.Nincr,1); % Substrate spring energy
    Energy.Es=0;
    Energy.Ev=0;
    Energy.Ef=0;
    Energy.Eb=0;
    Energy.EB=0;
    Energy.Ec=0;
    Energy.Ei=0;
    Energy.Esub=0;
elseif nargin>2
    Energies.S(incr)=Energy.Es;
    Energies.V(incr)=Energy.Ev;
    Energies.B(incr)=Energy.EB;
    if Set.Bending
        Energies.b(incr)=Energy.Eb;
    end
    Energies.F(incr)=Energy.Ef;
    if Set.Contractility && (Set.cPurseString > 0 || Set.cLateralCables > 0)
        Energies.C(incr)=Energy.Ec;
    end
    % Changed by Adria. Controversial. Currently undefined, have to see why
    % is so.
%     if Set.Substrate
%         Energies.Sub(incr) = Energy.Esub;
%     end
end
end