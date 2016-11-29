function [ output_args ] = constellation(points)

bits = de2bi(0:length(points)-1);
bits = fliplr(bits);
%scatterplot(points);

scatter(real(points), imag(points),'+r');
axis([-1.15 1.15 -1.15 1.15]);
for j = 1:length(points)
b = num2str(bits(j,:));
c = cellstr(b);
d = strrep(c, ' ', '');
text(real(points(j))-.1, imag(points(j))+.06, d, 'FontSize',9);
end

end

