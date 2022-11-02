function [Ynew, Tnew] = YFlip6N(Ys, Ts, XsToDisconnect, Geo)
%YFLIP6N Summary of this function goes here
%   Detailed explanation goes here
Tnew = [];
Ynew = [];

Xs = unique(Ts);
Xs_g = Xs(ismember(Xs, Geo.XgID));
Xs_c = Xs(~ismember(Xs, Geo.XgID));

[Xs_gConnectedNodes, Xs_gUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_g);
[Xs_cConnectedNodes, Xs_cUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_c);

end

