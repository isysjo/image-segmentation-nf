import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Set;
        
import javax.imageio.ImageIO;


public class LP3 {
	public static double getColorDistance(Color c1, Color c2) {
		int r1 = c1.getRed();
		int g1 = c1.getGreen();
		int b1 = c1.getBlue();
		
		int r2 = c2.getRed();
		int g2 = c2.getGreen();
		int b2 = c2.getBlue();
		
		return Math.sqrt((r1-r2)*(r1-r2) + (g1-g2)*(g1-g2) + (b1-b2)*(b1-b2));
	}
	
	public static Color averageColor(List<Color> colors) {
		int n = colors.size();
		int avgRed = 0;
		int avgGreen = 0;
		int avgBlue = 0;
		for(Color c : colors) {
			avgRed+=c.getRed();
			avgGreen+=c.getGreen();
			avgBlue+=c.getBlue();
		}
		avgRed/=n;
		avgGreen/=n;
		avgBlue/=n;
		return new Color(avgRed,avgGreen,avgBlue);
	}
	
	public static int n;
	
	public static void main(String[] args) {
		BufferedImage img = null;
              
		try {
		    img = ImageIO.read(new File("C:\\Users\\JTwal\\Documents\\NetBeansProjects\\LP3\\sagegrouse.jpg"));
		} catch (IOException e) {
			System.out.println("Error parsing image. " + e);
			return;
		}
		
		//Parse the image into a 2d array of colors
		int height = img.getHeight();
		int width = img.getWidth();

		n = height*width;
		
		System.out.println("Number of pixels: " + n + ", Image height: " + height + ", Image width: " + width);
		
		Color[][] colors = new Color[width][height];
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				int index = i*height + j;
				colors[i][j] = new Color(img.getRGB(i, j));
			}
		}

		// Begin two-means clustering of image pixels
		Color c1 = Color.WHITE;
		Color c2 = Color.BLACK;
		// Repeatedly set c1 and c2 to the average of the pixels that are closest in terms of "color"
                
		for(int t = 0; t < 20; t++) {
			List<Color> closeToC1 = new ArrayList<Color>(); 
			List<Color> closeToC2 = new ArrayList<Color>(); 
			for(int i = 0; i < width; i++) {
				for(int j = 0; j < height; j++) {
					double distToC1 = getColorDistance(c1,colors[i][j]);
					double distToC2 = getColorDistance(c2,colors[i][j]);
					
					if(distToC1 > distToC2) {
						closeToC1.add(colors[i][j]);
					} else {
						closeToC2.add(colors[i][j]);
					}
				}
			}
			
			c1 = averageColor(closeToC1);
			c2 = averageColor(closeToC2);
		}
		
		//TODO: Use flow based algorithm to improve the two-means image segmentation.
		// create graph (G,E) using pixels as V and neighboring pixels as E
                // calculate likelihoods a_v  and b_v for each pixel
                // vertices stored as their coordinates pairs
                // edges stored as pairs of vertices
                HashMap<Pair<Integer>,List<Integer>> graph = new HashMap();
                HashMap<Integer,Pair<Integer>> vertices = new HashMap();
                int a = 9;
                int b = 8;
                buildPixelGraph(colors,graph,c1,c2,a,b,vertices);
               
                // run flow algorithm
                double[][] finalRGraph = mincut(0,n+1,graph,vertices,a,b);
                
                // assign pixels to S and T sets
                HashSet<Integer> finalCluster1 = new HashSet();
                bfs(0,n+1,finalCluster1,finalRGraph,graph,vertices, new HashMap<Integer,Integer>());
               
		//Use c1 and c2 to decide which cluster each pixel belongs to
		BufferedImage outImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

                for(int vert : vertices.keySet()) {
                    if(vert!= 0 && vert!= n+1) {
                    if(finalCluster1.contains(vert))
                       outImage.setRGB(vertices.get(vert).first, vertices.get(vert).second, c1.getRGB());
                    else outImage.setRGB(vertices.get(vert).first, vertices.get(vert).second, c2.getRGB());}
                   
                }

		//Write the output image
		try {
			File outFile = new File("out_sage.png");
			ImageIO.write(outImage, "png", outFile);
		} catch (IOException e) {
			System.out.println("Error when writing output image. " + e);
		}
	}
        
        public static double[][]  mincut(int s, int t, HashMap<Pair<Integer>,List<Integer>> graph, 
                HashMap<Integer,Pair<Integer>> vertices,
                int a, int b)  {
            int u,v;
            double[][] rGraph = new double[graph.size()][graph.size()];  

            for (int i = 0; i < graph.size(); i++) { 
            for (int j = 0; j < graph.size(); j++) { 
                if(graph.get(vertices.get(i)).contains(j)) {
                      if(i == 0)
                        rGraph[i][j] = a;
                      else if(j == graph.size()-1)
                          rGraph[i][j] = b;
                      else rGraph[i][j] = 1;
                }       
                else rGraph[i][j] = 0;
               
            }
        }
            
            HashMap<Integer,Integer> pathSet = new HashMap();
            HashSet<Integer> visit = new HashSet();
        
            while (bfs(s, t,visit,rGraph,graph,vertices,pathSet)) { 
       
            // Find minimum residual capacity of the edges  
            // along the path filled by BFS.
            double pathFlow = Integer.MAX_VALUE;          
            for (v = t; v != s; v = pathSet.get(v)) {
                u = pathSet.get(v); 
                pathFlow = Math.min(pathFlow, rGraph[u][v]); 
                
            } 
            

            // update residual capacities of the edges and  
            // reverse edges along the path 
           
            for (v = t; v != s; v = pathSet.get(v)) {
                u = pathSet.get(v);  
                rGraph[u][v] = rGraph[u][v] - pathFlow; 
                rGraph[v][u] = rGraph[v][u] + pathFlow; 
            } 
            
            
            pathSet.clear();
            visit.clear();
        } 
        return rGraph;
        
        }
        
        public static boolean bfs(Integer start, Integer goal, Collection visited, 
                double[][] rGraph, HashMap<Pair<Integer>,List<Integer>> graph, 
                HashMap<Integer,Pair<Integer>> vertices,
                HashMap<Integer,Integer> pathSet) {
           
            if(start.equals(goal)) {pathSet.put(start,null); return true;}
            Queue<Integer> q = new LinkedList<Integer>();
            q.add(start);
            while(!q.isEmpty()) {
                
                int curr = q.remove();
                
                for(Object neighbor : graph.get(vertices.get(curr))) {
                    if(!visited.contains((Integer)neighbor) && 
                            rGraph[curr][(Integer)neighbor] > 0) {
                        visited.add((Integer)neighbor);
                        q.add((Integer)neighbor);
                        pathSet.put((Integer)neighbor, curr);
                        if(neighbor.equals(goal)) return true;
                    }
                }
              
            }
            
           return false; 
        }
        
        public static void buildPixelGraph(Color[][] image, HashMap<Pair<Integer>, List<Integer>> graph, 
                Color fg, Color bg, double a, double b, HashMap<Integer,Pair<Integer>> vertices) {
            Pair<Integer> source = new Pair(-1,-1);
            Pair<Integer> sink = new Pair(image.length,image[0].length);
            graph.put(source, new ArrayList<Integer>());
            graph.put(sink, new ArrayList<Integer>());
            vertices.put(0,source);
            vertices.put(image.length*image[0].length+1,sink);
        
            int index = 1;
            for(int i = 0; i < image.length; i++) {
                for(int j = 0; j < image[0].length; j++) {
                    Pair<Integer> v = new Pair(i,j);
                    graph.put(v, new ArrayList<Integer>());
                    vertices.put(index, v);
                    
                    //connect to source or sink
                    if(getColorDistance(fg,image[i][j]) > getColorDistance(bg,image[i][j])) {
                        List c1 = (ArrayList)graph.get(source);
                        c1.add(index);    
                     
                                
                    } else {
                        List c2 = (ArrayList)graph.get(v);
                        c2.add(image.length*image[0].length+1);
                     
                    }
                    
                    //find adjacent pixels
                    Pair<Integer> u;
                    List adj = graph.get(v);
                    if(j + 1 < image[0].length) {
                        u = new Pair<Integer>(i,j+1);
                   
                        adj.add(u.first*image[0].length+u.second+1);
                        graph.put(v, adj);
                
                    }
                    
                    if(j - 1 >= 0) {
                        u = new Pair<Integer>(i,j-1);
           
                        adj.add(u.first*image[0].length+u.second+1);
                        graph.put(v, adj);
                 
                    }
                    
                    if(i + 1 < image.length) {
                        u = new Pair<Integer>(i+1,j);
              
                         adj.add(u.first*image[0].length+u.second+1);
                        graph.put(v, adj);
                     
                    }
                    
                    if(i - 1 >= 0) {
                        u = new Pair<Integer>(i-1,j);
                     
                        adj.add(u.first*image[0].length+u.second+1);
                        graph.put(v, adj);
                      
                    }
                    
                    if(j + 1 < image[0].length && i - 1 >= 0) {
                        u = new Pair<Integer>(i-1,j+1);
             
                       adj.add(u.first*image[0].length+u.second+1);
                        graph.put(v, adj);
                      
                    }
                    
                    if(j + 1 < image[0].length && i +1 < image.length) {
                        u = new Pair<Integer>(i+1,j+1);
                        adj.add(u.first*image[0].length+u.second+1);
                        graph.put(v, adj);
        
                    }
                    
                    if(j - 1 >= 0 && i - 1 >= 0) {
                        u = new Pair<Integer>(i-1,j-1);
                       adj.add(u.first*image[0].length+u.second+1);
                        graph.put(v, adj);
                      
                    }
                    
                    if(j - 1 >= 0 && i + 1 < image.length) {
                       u = new Pair<Integer>(i+1,j-1);
                        adj.add(u.first*image[0].length+u.second+1);
                        graph.put(v, adj);
                    }
                    
                    index++;    
                }
            }
            
         
        }
        
        public static class Pair<T> {
            T first;
            T second;
            public Pair(T f, T s) {
                first = f;
                second = s;
            }
            
            public boolean equals(Pair<T> p2) {
                return first.equals(p2.first) && second.equals(p2.second);
            }
            
            public String toString() {
                return "first: " + first + " second: " + second;
            }
}


}
