package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class DefaultTeam {

	private int compute_quadrant(Point p, Point maison, ArrayList<Integer> alphas) {
		double m = (float)(p.getY() - maison.getY()) / (p.getX() - maison.getX());
		
		double angle = Math.toDegrees(Math.atan(m));
		
		if (angle > 0.0) {
			if (p.getX() < maison.getX()) { // Third quadrant
				angle = 180 + angle;
			}
		} else {
			if (p.getX() < maison.getX()) { // Second quadrant
				angle = 180 + angle;
			} else { // fourth quadrant
				angle = 360 + angle;
			}
		}
		
		for (int i=0; i < alphas.size(); i++) {
			if (i == alphas.size()-1) {
				return i;
			}

			if (angle >= alphas.get(i) && angle <= alphas.get(i+1)) {
				return i;
			}
		}
		
		return 0;
	}


	private double score(ArrayList<Point> points) {
    	double score = 0;
    	for (int i=1; i<points.size(); i++) {
    		score += points.get(i).distance(points.get(i-1));
		}
		
		score += points.get(0).distance(points.get(points.size()-1));
  	
    	return score;
    }
	
	
    private ArrayList<Point> improve(ArrayList<Point> points) {    	
    	for (int i=0;i<points.size();i++){
            for (int j=i+2;j<points.size() ;j++){
                double a=points.get(i).distance(points.get((i+1)%points.size()));
                double b=points.get(j%points.size()).distance(points.get((j+1)%points.size()));
                double c=points.get(i).distance(points.get(j%points.size()));
                double d=points.get((i+1)%points.size()).distance(points.get((j+1)%points.size()));
                if (a+b>c+d) {
                    ArrayList<Point> p=new ArrayList<Point>();
                    for (int k=0;k<=i;k++) p.add(points.get(k));
                    for (int k=j;k>i;k--) p.add(points.get(k));
                    for (int k=j+1;k<points.size();k++) p.add(points.get(k));
                    return p;
                }
            }
        }
        return points;
	}


	public ArrayList<Point> calculTSP(ArrayList<Point> points, Point startpoint) {
		if (points.size()<4) {
            return points;
        }

        ArrayList<Point> return_points = new ArrayList<Point>();
        ArrayList<Point> rest = (ArrayList<Point>) points.clone();

        return_points.add(startpoint);

        while(!rest.isEmpty()) {
            Point lastVisited = return_points.get(return_points.size()-1);
            Point nearest = rest.get(0);

            for(Point r:rest) {
                if(lastVisited.distance(r) < lastVisited.distance(nearest))
                    nearest=r;
            }
            return_points.add(nearest);
            rest.remove(nearest);
        }
        
        // Now let's improve the solution by removing the overlapping paths
        // Local optimization
        double current_score = score(return_points);
        double old_score = current_score +1;
        
        while (old_score > current_score) {
        	return_points = improve(return_points);
        	
        	old_score = current_score;
        	current_score = score(return_points);
        }

        return return_points;
	}

	
	private ArrayList<Point> reduce_budget(ArrayList<Point> points, Point maison, double budget) {
		double currentCost = score(points);
		Point pointToRemove;

		TreeMap<Double, ArrayList<Point>> costDecreaseHeap = new TreeMap<Double, ArrayList<Point>>();
		TreeMap<Integer, Double> reverseCostDecreaseHeap = new TreeMap<Integer, Double>();

		for (Point p: points) {
			int index = points.indexOf(p);

			if (index == points.indexOf(maison)) continue;

			Point pointBefore = points.get((index-1)%points.size());
			Point pointAfter = points.get((index+1)%points.size());

			double distance_a = p.distance(pointBefore);
			double distance_b = p.distance(pointAfter);
			double distance_x = pointBefore.distance(pointAfter);

			double costDecrease = ((distance_a + distance_b) - distance_x);

			if (costDecreaseHeap.containsKey(costDecrease)) costDecreaseHeap.get(costDecrease).add(p);
			else costDecreaseHeap.put(costDecrease, new ArrayList<Point>(Arrays.asList(p)));

			reverseCostDecreaseHeap.put(index,costDecrease);
		}

		while (currentCost > budget) {
			double mostCostDecrease = costDecreaseHeap.lastKey();
			if (costDecreaseHeap.get(mostCostDecrease).size() > 1) {
				pointToRemove = costDecreaseHeap.get(mostCostDecrease).remove(0);
			}
			else {
				pointToRemove = costDecreaseHeap.remove(mostCostDecrease).get(0);
			}

			int indexOfPointToRemove = points.indexOf(pointToRemove);
			points.remove(indexOfPointToRemove);

			currentCost -= mostCostDecrease;

			double old_score = currentCost + 1;
			while (old_score > currentCost) {
				points = improve(points);
				
				old_score = currentCost;
				currentCost = score(points);
			}			
			
			Point pointBefore = points.get((indexOfPointToRemove-1)%points.size());
			Point pointAfter = points.get((indexOfPointToRemove+1)%points.size());

			double pointBeforeCostDecrease = reverseCostDecreaseHeap.get(indexOfPointToRemove-1);
			double pointAfterCostDecrease = reverseCostDecreaseHeap.get(indexOfPointToRemove+1);

			if (costDecreaseHeap.get(pointBeforeCostDecrease) != null) {
				if (costDecreaseHeap.get(pointBeforeCostDecrease).size() > 1) {
					costDecreaseHeap.get(pointBeforeCostDecrease).remove(costDecreaseHeap.get(pointBeforeCostDecrease).indexOf(pointBefore));
				}
				else {
					costDecreaseHeap.remove(pointBeforeCostDecrease).get(0);
				}
			}

			if (costDecreaseHeap.get(pointAfterCostDecrease) != null) {
				if (costDecreaseHeap.get(pointAfterCostDecrease).size() > 1) {
					costDecreaseHeap.get(pointAfterCostDecrease).remove(costDecreaseHeap.get(pointAfterCostDecrease).indexOf(pointAfter));
				}
				else {
					costDecreaseHeap.remove(pointAfterCostDecrease).get(0);
				}
			}

			reverseCostDecreaseHeap.remove(indexOfPointToRemove-1);
			reverseCostDecreaseHeap.remove(indexOfPointToRemove+1);

			Point pointBeforePointBefore = points.get((indexOfPointToRemove-1-1)%points.size());
			Point pointAfterPointBefore = points.get((indexOfPointToRemove-1+1)%points.size());

			double distance_a = pointBefore.distance(pointBeforePointBefore);
			double distance_b = pointBefore.distance(pointAfterPointBefore);
			double distance_x = pointBeforePointBefore.distance(pointAfterPointBefore);

			double pointBeforeNewCostDecrease = ((distance_a + distance_b) - distance_x);

			if (costDecreaseHeap.containsKey(pointBeforeNewCostDecrease)) costDecreaseHeap.get(pointBeforeNewCostDecrease).add(pointBefore);
			else costDecreaseHeap.put(pointBeforeNewCostDecrease, new ArrayList<Point>(Arrays.asList(pointBefore)));

			reverseCostDecreaseHeap.put(indexOfPointToRemove-1,pointBeforeNewCostDecrease);

			Point pointBeforePointAfter = points.get((indexOfPointToRemove+1-1)%points.size());
			Point pointAfterPointAfter = points.get((indexOfPointToRemove+1+1)%points.size());

			distance_a = pointAfter.distance(pointBeforePointAfter);
			distance_b = pointAfter.distance(pointAfterPointAfter);
			distance_x = pointBeforePointAfter.distance(pointAfterPointAfter);

			double pointAfterNewCostDecrease = ((distance_a + distance_b) - distance_x);

			if (costDecreaseHeap.containsKey(pointAfterNewCostDecrease)) costDecreaseHeap.get(pointAfterNewCostDecrease).add(pointAfter);
			else costDecreaseHeap.put(pointAfterNewCostDecrease, new ArrayList<Point>(Arrays.asList(pointAfter)));

			reverseCostDecreaseHeap.put(indexOfPointToRemove+1,pointAfterNewCostDecrease);
		}
		
		return points;
	}


  	public ArrayList<ArrayList<Point>> calculCinqVoyageursAvecBudget(Point maison, ArrayList<Point> points) {
		  
		double budget = 1664.0;

		TreeMap<Integer, ArrayList<ArrayList<Point>>> solutions = new TreeMap<Integer, ArrayList<ArrayList<Point>>>();

		for (int alpha = 0; alpha < 72; alpha++) {
    
			ArrayList<Point> alice = new ArrayList<Point>(); alice.add(maison);
			ArrayList<Point> bob = new ArrayList<Point>(); bob.add(maison);
			ArrayList<Point> cindy = new ArrayList<Point>(); cindy.add(maison);
			ArrayList<Point> dave = new ArrayList<Point>(); dave.add(maison);
			ArrayList<Point> eddy = new ArrayList<Point>(); eddy.add(maison);
			
			ArrayList<Point> unvisited = (ArrayList<Point>) points.clone();
			unvisited.remove(maison);
			
			
			// // Removing the points too far away
			// for (Point p : unvisited) {
			// 	if (p.distance(maison) >= budget/2) {
			// 		unvisited.remove(p);
			// 	}
			// }
			

			// Divide the space in 5 parts
			ArrayList<Integer> alphas = new ArrayList<Integer>();

			for (int i = 0; i < 5; i++) {
				alphas.add((i*72)+alpha);
			}

			Map<Integer, ArrayList<Point>> div_points = new HashMap<>();
			
			for (Point p : unvisited) {
				int quadrant = compute_quadrant(p, maison, alphas);
				ArrayList<Point> res = div_points.putIfAbsent(quadrant, new ArrayList<Point>(Arrays.asList(p)));
				if (res != null) {
					res.add(p);
					div_points.put(quadrant, res);
				}
			}

			
			// For each person find the total traversed path
			alice = calculTSP(div_points.get(0), maison);
			bob = calculTSP(div_points.get(1), maison);
			cindy = calculTSP(div_points.get(2), maison);
			dave = calculTSP(div_points.get(3), maison);
			eddy = calculTSP(div_points.get(4), maison);
			
			
			// // For each person, reduce it under the budget
			alice = reduce_budget(alice, maison, budget);
			bob = reduce_budget(bob, maison, budget);
			cindy = reduce_budget(cindy, maison, budget);
			dave = reduce_budget(dave, maison, budget);
			eddy = reduce_budget(eddy, maison, budget);
			
			ArrayList<ArrayList<Point>> result = new ArrayList<ArrayList<Point>>();
			result.add(alice);
			result.add(bob);
			result.add(cindy);
			result.add(dave);
			result.add(eddy);

			int resultSize = 0;

			for (int i = 0; i < result.size(); i++) resultSize += result.get(i).size();

			solutions.put(resultSize, result);
		}
	
		return solutions.get(solutions.lastKey());
  	}
}