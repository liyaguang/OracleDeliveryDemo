package edu.usc.infolab.main;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import javax.management.Query;

import edu.usc.infolab.lib.GeoPoint;
import edu.usc.infolab.lib.QueryBuilder;

public class DataRetrievalTest {
	public static void main(String[] argv) throws IOException {
		// List<GeoPoint> points = getNodes();
		// genDistFile(points);
	}

	public static void genDistFile(List<GeoPoint> points)
			throws IllegalStateException, IOException {
		String logFileName = "data/history.txt";
		String outFileName = "data/dist.txt";
		String day = "Tuesday";
		int time = 30;
		int size = points.size();
		PrintWriter pr = new PrintWriter(outFileName, "UTF-8");
		PrintWriter logWriter = new PrintWriter(logFileName, "UTF-8");
		for (int s = 0; s < size; ++s) {
			for (int e = 0; e < size; ++e) {
				String travelTime = "0";
				if (s != e) {
					String result = QueryBuilder.getResult(points.get(s),
							points.get(e), time, day);
					logWriter.println(String.format("n%d->n%d\t%s", s, e,
							result));
					String[] fs = result.split("@");
					if (fs.length > 1) {
						travelTime = fs[0]
								.substring(fs[0].lastIndexOf(';') + 1);
						travelTime = travelTime.substring(0,
								travelTime.indexOf('-'));

					}
				}
				pr.write(String.format("%.3f", Double.parseDouble(travelTime)));
				if (e != size - 1) {
					pr.write("\t");
				}
			}
			pr.println();
		}
		logWriter.close();
		pr.close();
	}

	public static double[][] getDistArray() throws IOException {
		// read distance array
		String filename = "data/dist.txt";
		List<String> lines = Files.readAllLines(Paths.get(filename),
				StandardCharsets.UTF_8);
		int count = lines.size();
		double[][] array = new double[count][count];
		for (int i = 0; i < count; ++i) {
			String[] fs = lines.get(i).split("\t");
			for (int j = 0; j < fs.length; ++j) {
				array[i][j] = Double.parseDouble(fs[j]);
			}
		}
		return array;
	}

	public static List<GeoPoint> getNodes() throws IOException {
		String nodeFileName = "data/Nodes.csv";
		String sampledNodeFileName = "data/sample_nodes.csv";
		PrintWriter pw = new PrintWriter(sampledNodeFileName, "UTF-8");
		List<String> lines = Files.readAllLines(Paths.get(nodeFileName),
				StandardCharsets.UTF_8);
		// get sample points
		List<GeoPoint> points = new ArrayList<GeoPoint>();
		int stepSize = 200, count = 10;
		assert lines.size() > stepSize * count;
		for (int i = 0; i < count; ++i) {
			String line = lines.get(i * stepSize);
			String[] fields = line.split(",");
			double lat = Double.parseDouble(fields[1]);
			double lng = Double.parseDouble(fields[2]);
			GeoPoint point = new GeoPoint(lat, lng);
			points.add(point);
			pw.println(String.format("n%d,%s", i, point.toString()));
		}
		pw.close();
		return points;
	}
}
