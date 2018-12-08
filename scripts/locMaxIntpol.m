function [locMax] = locMaxIntpol(Z, cut)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Set Default values
if(nargin==1)
    cut = 1;
end

% Get size of the field
sZ  = size(Z);

Zp = interpn(Z, 'makima');

% Find local maxima of the field Z and the indices
Imax  = imregionalmax(Z);
Imaxp  = imregionalmax(Zp);

Zmax  = Z(find(Imax==1));
Zpmax = Zp(find(Imaxp==1));

[Zmax'; Zpmax']

[Ix,Iy,Iz] = ind2sub(sZ, find(Imax==1));
        
% remove the maxima on the boundary
% z-boundaries
z_valid = ~(Iz == sZ(3) | Iz == 1 ) & ~(Iy == sZ(2) | Iy == 1 ) & ~(Ix == sZ(2) | Ix == 1 );
Ix = Ix(z_valid);
Iy = Iy(z_valid);
Iz = Iz(z_valid);

Zneighbours = zeros([1 27]);

% Fit quadratic to the data
modelfun = @(b,x) b(1) + sum(x * reshape(b(2:10), [3,3]) .* x,2);

modelfun = @(b,x) b(1) + sum( (x-[b(2) b(3) b(4)]) * [[b(5) b(8)  b(9)];...
                                                      [b(8) b(6)  b(10)];...
                                                      [b(9) b(10) b(7)]] .* (x-[b(2) b(3) b(4)]),2);
                                                  
modelfun = @(b,x) b(1) + sum( (x-[b(2) b(3) b(4)]) * diag(b(5:7)) .* (x-[b(2) b(3) b(4)]),2);

modelfun = @(b,x) b(1) + sum( (x-x(1,:)) * [[b(2) b(5)  b(6)];...
                                            [b(5) b(3)  b(7)];...
                                            [b(6) b(7)  b(4)]] .* (x-x(1,:)),2);
                                        
modelfun = @(b,x) b(1) + sum( (x-x(1,:)) * diag(b(2:4)) .* (x-x(1,:)),2);
                                        
                                        
locMax1 = zeros([1, length(Ix)]);
locMax  = zeros([1, length(Ix)]);
% Fit the maxima using a quadratic to a 26 neighbourhood of the discrete
% maxima
for i = 1:length(Ix)
    %tmp1 = sub2ind(sZ, Ix(i), Iy(i), Iz(i))
    
 %   [Iadj , Radj, Nfound ] = neighbourND( tmp1, sZ, [1 1 1] );
    
    % indices of discrete maximum and neighbourhood
    relevant_indices = [
                    [Ix(i)   Iy(i)   Iz(i)  ];
                    [Ix(i)   Iy(i)-1 Iz(i)  ];
                    [Ix(i)   Iy(i)+1 Iz(i)  ];
                    [Ix(i)   Iy(i)   Iz(i)-1];
                    [Ix(i)   Iy(i)   Iz(i)+1];
                    [Ix(i)   Iy(i)-1 Iz(i)+1];
                    [Ix(i)   Iy(i)+1 Iz(i)-1];
                    [Ix(i)   Iy(i)-1 Iz(i)-1];
                    [Ix(i)   Iy(i)+1 Iz(i)+1];
                    
        
                    [Ix(i)-1 Iy(i)   Iz(i)  ];
                    [Ix(i)-1 Iy(i)-1 Iz(i)  ];
                    [Ix(i)-1 Iy(i)+1 Iz(i)  ];
                    [Ix(i)-1 Iy(i)   Iz(i)-1];
                    [Ix(i)-1 Iy(i)   Iz(i)+1];
                    [Ix(i)-1 Iy(i)-1 Iz(i)+1];
                    [Ix(i)-1 Iy(i)+1 Iz(i)-1];
                    [Ix(i)-1 Iy(i)-1 Iz(i)-1];
                    [Ix(i)-1 Iy(i)+1 Iz(i)+1];
                    
                    [Ix(i)+1 Iy(i)   Iz(i)  ];
                    [Ix(i)+1 Iy(i)-1 Iz(i)  ];
                    [Ix(i)+1 Iy(i)+1 Iz(i)  ];
                    [Ix(i)+1 Iy(i)   Iz(i)-1];
                    [Ix(i)+1 Iy(i)   Iz(i)+1];
                    [Ix(i)+1 Iy(i)-1 Iz(i)+1];
                    [Ix(i)+1 Iy(i)+1 Iz(i)-1];
                    [Ix(i)+1 Iy(i)-1 Iz(i)-1];
                    [Ix(i)+1 Iy(i)+1 Iz(i)+1];
                    ];
                for j = 1:27
                    Zneighbours(j) = Z( sub2ind(sZ,relevant_indices(j,1), ...
                                        relevant_indices(j,2), relevant_indices(j,3)));
                end
                
                Zp = interpn(Z, 'cubic');
                
     %           Zneighbours2 = Z( Iadj );
                mdl = fitnlm( relevant_indices, Zneighbours', modelfun,...
                              [ Zneighbours(1) ones([1,6])] );
                locMax(i)  = mdl.Coefficients.Estimate(1);
                locMax1(i) = Zneighbours(1);
                
      %          [Iadj , Radj, Nfound ] = neighbourND( index, sizeA, res )
end

[locMax; locMax1]

Indices_max = find(Imax==1);

[Ix,Iy,Iz] = ind2sub(smask_valid, Indices_max);

if( minima==1 )
    Imin = imregionalmin(Z);
end

locMax = 
end

