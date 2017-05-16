args <- commandArgs(trailingOnly=T)

csv.file = 'res.csv'
pdf.file = 'mpc.pdf'

hlines = c()

if (length(args) >= 1) {
    csv.file = args[1]
    pdf.file = gsub('.csv', '.pdf', csv.file)
    if (length(args) > 1)
        hlines = args[-1]
}

# load csv output
d <- read.csv(csv.file)
# sampling time
ts = 0.1
# add time to dataset
d$time = 0:(nrow(d)-1) * ts

if (!interactive()) {
    pdf(pdf.file, width=16/2.8, height=9/2.8)
    par(mar=c(3.2, 3.2, 0.3, 0.3))
}

# create a plot
plot.new()
plot.window(xlim=c(0, max(d$time)), ylim=c(-1.6, 1.6), xaxs="i", yaxs="i")

lines(d$time, d$v, lwd=2, col=1, lty=1)
abline(h=1, lwd=2, col=2, lty=2)
lines(d$time, d$a, lwd=2, col=3, lty=3)
lines(d$time, d$u, lwd=2, col=4, lty=4)
lines(d$time, d$du, lwd=2, col=5, lty=5)
if (length(hlines) > 0) {
    for (h in hlines)
        abline(h=h, lwd=2, col=gray(0.5), lty=2)
}

legend("bottomright",
       c('speed (m/s)', 'reference (m/s)', 'acceleration (m/s/s)', 'control (m/s/s)', 'control derivative (m/s/s/s)'),
       ncol=1,
       inset=c(0, 0),
       lty=1:5,
       col=1:5,
       lwd=2,
       box.lwd=0,
       pch=NA,
       pt.bg='white',
       cex=.8,
       seg.len=3)

title(xlab='time (s)')
axis(1, las=1)
axis(2, las=1)
box()

if (!interactive())
    dev.off()
