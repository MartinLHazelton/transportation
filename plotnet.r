#' Network plotting function
#'
#' This function plots traffic networks, where the relative locations of nodes are meaningful (unlike conceptual networks).
#' @param nodes.x Vector of x-coordinates for node locations
#' @param nodes.y Vector of y-coordinates for node locations
#' @param link.node.matrix Two column matrix specifying nodes of origin and destination for each link
#' @param radius Radius of circles marking node locations. Defaults to 0.5.
#' @param link.sep Angular separation (as a fraction of pi radians) of entry/exit positions of bi-directional links. Defaults to 1/12.
#' @param link.label.offset Controls position of link labels relative to link. Detauls to 0.4.
#' @param LABEL.CEX Character expansion factor for link labels. Defaults to 1.
#' @param ARROW.LENGTH Controls length of arrow heads. Defaults to 0.4.
#' @param ARR.TYPE Controls appearance of arrow heads in Arrows function. Defaults to "triangle".
#' @param node.label Label nodes? Defaults to TRUE.
#' @param LINK.LABEL Label links? Defaults to TRUE.
#' @param RightAbove Further control for link label positions. Defaults to FALSE, which means that labels appear to left looking in direction of arrows. If TRUE, labels appear to the left/above links.
#' @param link.lty Line type for links: single value, or vector (lty per link). Defaults to 1.
#' @param link.col Line type for links: single value, or vector (lty per link). Defaults to 1 (black).
#' @param node.lwd Line width for nodes: single value, or vector (lty per nodes). Defaults to 2. 
#' @return Produces a plot of the network
#' @keywords network plotting
#' @export
#' @examples
#' nodes.x <- c(2,7,12,2,7,12)
#' nodes.y <- c(7,7,7,2,2,2)
#' link.node.matrix <- rbind(c(1,2),c(1,4),c(2,5),c(3,2),c(3,6),c(4,5),c(6,5))
#' plotnet(nodes.x,nodes.y,link.node.matrix,ARROW.LENGTH=0.3,link.sep=0,radius=0.7,RightAbove=T)

plotnet <- function(nodes.x,nodes.y,link.node.matrix,radius=1/2,link.sep=1/12,link.label.offset=0.4,LABEL.CEX=1,ARROW.LENGTH=0.4,ARR.TYPE="triangle",node.label=TRUE,LINK.LABEL=TRUE,RightAbove=F,link.lty=1,link.col=1,node.lwd=2){
  require(shape)
  n.nodes <- length(nodes.x)
  n.links <- nrow(link.node.matrix)
  if (length(link.sep) < n.links) link.sep <- rep(link.sep,length=n.links)
  if (length(link.lty) < n.links) link.lty <- rep(link.lty,length=n.links)
  if (length(link.col) < n.links) link.col <- rep(link.col,length=n.links)
  if (length(node.lwd) < n.nodes) node.lwd <- rep(node.lwd,length=n.nodes)
  plot(0,0,xlim=c(min(nodes.x)-radius,max(nodes.x)+radius),ylim=c(min(nodes.y)-radius,max(nodes.y)+radius),type="n",axes=F,xlab="",ylab="")
  for (i in 1:n.nodes){
	symbols(nodes.x[i],nodes.y[i],circles=radius,add=T,inches=FALSE,lwd=node.lwd[i]) 
  }
  if (node.label) text(nodes.x,nodes.y,as.character(1:n.nodes),cex=LABEL.CEX,font=2)
  for (i in 1:n.links){
    from.node <- link.node.matrix[i,1]
    to.node <- link.node.matrix[i,2]
    from.x <- nodes.x[from.node]
    from.y <- nodes.y[from.node]
    to.x <- nodes.x[to.node]
    to.y <- nodes.y[to.node]
    
    theta <- atan2((to.y-from.y),(to.x-from.x))
    
    theta1 <- theta + pi*link.sep[i]
    theta2 <- theta - pi*link.sep[i]
    
    y1 <- from.y + radius*1.1*sin(theta1)
    y2 <- to.y - radius*1.1*sin(theta2)
    x1 <- from.x + radius*1.1*cos(theta1)
    x2 <- to.x - radius*1.1*cos(theta2)
    
    Arrows(x1,y1,x2,y2,lwd=1,arr.adj=1,arr.length=ARROW.LENGTH,arr.type=ARR.TYPE,lty=link.lty[i],col=link.col[i])
    
    if (LINK.LABEL) {
      if (!RightAbove){
        x.pos <- (x1+x2)/2-link.label.offset*sin(theta1)
        y.pos <- (y1+y2)/2+link.label.offset*cos(theta2)
      }
      if (RightAbove){
        x.pos <- (x1+x2)/2+abs(link.label.offset*sin(theta1))
        y.pos <- (y1+y2)/2+abs(link.label.offset*cos(theta2))
      }     
      text(x.pos,y.pos,as.character(i),cex=LABEL.CEX,font=1)
    }
  }			
}