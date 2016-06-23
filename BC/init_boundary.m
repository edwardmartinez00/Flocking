function Bound=init_boundary(P)
    switch P.choiceBC
      case {1,2,3,4,5}
        theta = linspace(0,2*pi);
        BoundX = P.R*cos(theta);
        BoundY = P.R*sin(theta);
      case {6,7,8,9,10}
        BoundX = [linspace(-P.L,P.L), -P.L*ones(1,100), linspace(-P.L,P.L), P.L*ones(1,100)];
        BoundY = [-P.L*ones(1,100), linspace(-P.L,P.L), P.L*ones(1,100), linspace(-P.L,P.L)]; 
      otherwise
        BoundX = [];
        BoundY = [];
    end
    % finally
    Bound = [BoundX; BoundY];
end
