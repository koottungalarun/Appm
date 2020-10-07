
figure(1)
clf
hold on

phi = linspace(0,2*pi);
z = exp(-1i*phi);
plot(real(z), imag(z), 'k-')

showFaceNumber = true;
showVertexNumber = true;
for level = 0
    switch level
        case 0
            color = 'r';
        case 1
            color = 'g';
        case 2
            color = 'b';
    end
    
    filename = strcat('level', num2str(level), '-coords.dat');
    vc = importdata(filename);
    plot3(vc(:,1), vc(:,2), vc(:,3), 'o', 'Color', color)
    if showVertexNumber
        hold on
        text(vc(:,1), vc(:,2), vc(:,3), num2cell (1:size(vc,1)), 'Color', color)
        hold off
    end
    
    filename = strcat('level', num2str(level), '-f2v.dat');
    F = importdata(filename) + 1;
    patch('Faces', F, 'Vertices', vc, 'FaceColor', color, 'FaceAlpha', 0.1)
    shg
    axis equal
    
    if showFaceNumber
        % face centers
        nF = size(F,1);
        for i = 1 : nF
            fc(i,:) = mean(vc(F(i,:),:));
        end
        txt = num2cell(0:nF-1);
        text(fc(:,1), fc(:,2), fc(:,3), txt)
    end
end

shg
axis equal
hold off

