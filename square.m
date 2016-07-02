function [ boundary ] = square( L,C )
%UNTITLED4 Creates square to plot given side length and center
edgesX=[linspace(C(1)-L,C(1)+L),C(1)-L*ones(1,100),linspace(C(1)-L,C(1)+L),C(1)+L*ones(1,100)];
edgesY=[C(2)-L*ones(1,100),linspace(C(2)-L,C(2)+L),C(2)+L*ones(1,100),linspace(C(2)-L,C(2)+L)];
boundary=[edgesX;edgesY];

end

