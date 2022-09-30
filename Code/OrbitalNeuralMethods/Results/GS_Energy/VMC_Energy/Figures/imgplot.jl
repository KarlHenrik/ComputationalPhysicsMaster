import FileIO
import ColorVectorSpace
function imgplot(;file, extent, imgalpha, xticks, yticks)
    img = FileIO.load(file) .- ColorVectorSpace.RGBA(0, 0, 0, 1 - imgalpha)
    
    h, w = size(img)
    x0, x1, y0, y1 = extent
    dx = x1 - x0
    x_mid = (x0 + x1) / 2
    dy = y1 - y0
    
    plt.plot(img,
             xlim = (0, w),
             ylim = (0, h),
             yticks = h .- (yticks .- y0) .* h ./ dy,
             xticks = (xticks .+ dx ./ 2) .* w ./ dx,
             grid = true, gridalpha = 0.7,
             size = (600, 400),
             xformatter = x -> round((x / w - 0.5 + x_mid/dx) * dx, digits = 1),
             yformatter = y -> round((-y / h * dy + y0 + dy), digits = 10),
             )
    return x -> (x - x_mid + dx/2) * w / dx, y -> h - (y - y0) * h / dy
end