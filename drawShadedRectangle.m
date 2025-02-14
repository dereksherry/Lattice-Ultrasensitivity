function shapeHandle = drawShadedRectangle(axH,nameArray,valueArray,shadeType,fillColor,shadeColor,rotation)
if nargin < 6
    fillColor = 'none';
    shadeColor = [0,0,0];
end
if nargin < 4
    shadeType = 'none';
end
if nargin < 7
    rotation = 0;
end

shapeHandle = hgtransform;
r1 = rectangle(axH,nameArray,valueArray,'FaceColor',fillColor);
x0 = r1.Position(1);
y0 = r1.Position(2);
r1.Position = [-r1.Position(3)/2,-r1.Position(4)/2,r1.Position(3),r1.Position(4)];
r1.Parent = shapeHandle;

if ~strcmp('none',shadeType)
    hold(axH,'on')
    dotDensity = 7;
    r = r1.Curvature.*r1.Position(3:4)/2;
    counter = 0;
    dotCoordinates = [];
    
    for dx = 1/dotDensity:1/dotDensity:r1.Position(3)-1/dotDensity
        x = r1.Position(1) + dx;
        for dy = 1/dotDensity:1/dotDensity:r1.Position(4)-1/dotDensity
            y = r1.Position(2) + dy;
            if dx < r(1)
                xrel = r(1)-dx;
                if dy < r(2) % bottom left
                    yrel = r(2)-dy;
                    if (xrel/r(1))^2+(yrel/r(2))^2 < 1
                        counter = counter+1;
                        dotCoordinates(1,counter) = x;
                        dotCoordinates(2,counter) = y;
                    end
                elseif dy < r1.Position(4)-r(2) % middle left
                    counter = counter+1;
                    dotCoordinates(1,counter) = x;
                    dotCoordinates(2,counter) = y;
                else % top left
                    yrel = dy - (r1.Position(4) - r(2));
                    if (xrel/r(1))^2+(yrel/r(2))^2 < 1
                        counter = counter+1;
                        dotCoordinates(1,counter) = x;
                        dotCoordinates(2,counter) = y;
                    end
                end
            elseif dx <  r1.Position(3) - r(1) % middle, all height
                counter = counter+1;
                dotCoordinates(1,counter) = x;
                dotCoordinates(2,counter) = y;
            else
                xrel = dx - (r1.Position(3) - r(1));
                if dy < r(2) % bottom right
                    yrel = r(2)-dy;
                    if (xrel/r(1))^2+(yrel/r(2))^2 < 1
                        counter = counter+1;
                        dotCoordinates(1,counter) = x;
                        dotCoordinates(2,counter) = y;
                    end
                elseif dy < r1.Position(4)-r(2) % middle right
                    counter = counter+1;
                    dotCoordinates(1,counter) = x;
                    dotCoordinates(2,counter) = y;
                else % top right
                    yrel = dy - (r1.Position(4) - r(2));
                    if (xrel/r(1))^2+(yrel/r(2))^2 < 1
                        counter = counter+1;
                        dotCoordinates(1,counter) = x;
                        dotCoordinates(2,counter) = y;
                    end
                end
            end
        end
    end
    
    if strcmp('dot',shadeType) && ~isempty(dotCoordinates)
        shading = plot(dotCoordinates(1,:),dotCoordinates(2,:),'.','Color',shadeColor);
        shading.Parent = shapeHandle;
    end
    
    if strcmp('line',shadeType) && ~isempty(dotCoordinates)
        shading = plot(dotCoordinates(1,:),dotCoordinates(2,:),'Color',shadeColor);
        shading.Parent = shapeHandle;
    end
    
    counter = 0;
    startIndex = 1;
    if strcmp('vertical',shadeType) && ~isempty(dotCoordinates)
        for i = 2:size(dotCoordinates,2)
            if dotCoordinates(1,i) ~= dotCoordinates(1,i-1)
                counter = counter + 1;
                shading(counter) = plot(dotCoordinates(1,startIndex:i-1),dotCoordinates(2,startIndex:i-1),'Color',shadeColor);
                startIndex = i;
                shading(counter).Parent = shapeHandle;
            end
        end
        counter = counter + 1;
        i = i+1;
        shading(counter) = plot(dotCoordinates(1,startIndex:i-1),dotCoordinates(2,startIndex:i-1),'Color',shadeColor);
        shading(counter).Parent = shapeHandle;
    end
    
    
    if strcmp('horizontal',shadeType) && ~isempty(dotCoordinates)
        [~, newOrder] = sort(dotCoordinates(2,:));
        dotCoordinates = dotCoordinates(:,newOrder);
        for i = 2:size(dotCoordinates,2)
            if dotCoordinates(2,i) ~= dotCoordinates(2,i-1)
                counter = counter + 1;
                shading(counter) = plot(dotCoordinates(1,startIndex:i-1),dotCoordinates(2,startIndex:i-1),'Color',shadeColor);
                startIndex = i;
                shading(counter).Parent = shapeHandle;
            end
        end
        counter = counter + 1;
        i = i+1;
        shading(counter) = plot(dotCoordinates(1,startIndex:i-1),dotCoordinates(2,startIndex:i-1),'Color',shadeColor);
        shading(counter).Parent = shapeHandle;
    end
    
    if strcmp('slopeDown',shadeType) && ~isempty(dotCoordinates)
        [~, newOrder] = sort(dotCoordinates(2,:) + dotCoordinates(1,:));
        dotCoordinates = dotCoordinates(:,newOrder);
        for i = 2:size(dotCoordinates,2)
            if dotCoordinates(2,i) + dotCoordinates(1,i) ~= dotCoordinates(2,i-1) + dotCoordinates(1,i-1)
                counter = counter + 1;
                shading(counter) = plot(dotCoordinates(1,startIndex:i-1),dotCoordinates(2,startIndex:i-1),'Color',shadeColor,'linewidth',2);
                startIndex = i;
                shading(counter).Parent = shapeHandle;
            end
        end
        counter = counter + 1;
        i = i+1;
        shading(counter) = plot(dotCoordinates(1,startIndex:i-1),dotCoordinates(2,startIndex:i-1),'Color',shadeColor);
        shading(counter).Parent = shapeHandle;
    end
    
    if strcmp('slopeUp',shadeType) && ~isempty(dotCoordinates)
        [~, newOrder] = sort(abs(dotCoordinates(2,:) - dotCoordinates(1,:)));
        dotCoordinates = dotCoordinates(:,newOrder);
        for i = 2:size(dotCoordinates,2)
            if abs(dotCoordinates(2,i) - dotCoordinates(1,i)) ~= abs(dotCoordinates(2,i-1) - dotCoordinates(1,i-1))
                counter = counter + 1;
                shading(counter) = plot(dotCoordinates(1,startIndex:i-1),dotCoordinates(2,startIndex:i-1),'Color',shadeColor);
                startIndex = i;
                shading(counter).Parent = shapeHandle;
            end
        end
        counter = counter + 1;
        i = i+1;
        shading(counter) = plot(dotCoordinates(1,startIndex:i-1),dotCoordinates(2,startIndex:i-1),'Color',shadeColor);
        shading(counter).Parent = shapeHandle;
    end
end

R = makehgtform('zrotate',rotation);
T = makehgtform('translate',[x0,y0,0]);
shapeHandle.Matrix = shapeHandle.Matrix*T*R;
% shapeHandle.Matrix = shapeHandle.Matrix*t;
end