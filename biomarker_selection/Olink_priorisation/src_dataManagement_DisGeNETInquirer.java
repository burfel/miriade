/**
 * 
 */
package dataManagement;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.io.ObjectInputStream.GetField;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.HashMap;

import org.jdom2.Document;
import org.jdom2.JDOMException;

import webTools.HttpConnectionTools;

/**
 * @since 0.1.0 (15.07.2020)
 * @author Felicia Burtscher
 *
 */
public class DisGeNETInquirer {


	/**
	 * @since 0.1.0 (17.07.2020)
	 * @return
	 * @throws IOException 
	 * @throws JDOMException 
	 */
	public static int countUniprotHits(String uniprot) throws JDOMException, IOException {
		String queryURL = "https://www.disgenet.org/api/gda/gene/uniprot/" + uniprot;
		/*
		 * thanks to http://www.xyzws.com/Javafaq/how-to-use-httpurlconnection-post-data-to-web-server/139
		 */
		String queryParameters = null;
		try {
			queryParameters = "source=" + URLEncoder.encode("ALL", "UTF-8") +
					"&disease_class=" + URLEncoder.encode("C10", "UTF-8") +
					"&format=" + URLEncoder.encode("xml", "UTF-8");
		} catch (UnsupportedEncodingException e1) {
			// Auto-generated catch block
			e1.printStackTrace();
		}
		HashMap<String, String> requestMap = new HashMap<String, String>();
		requestMap.put("Content-Type", "application/xml");
		Document doc = 
				miscTools.XMLAuxTools.parseHttpRequestToXML(
						HttpConnectionTools.performGetRequest(queryURL, queryParameters, requestMap));
		return miscTools.XMLAuxTools.countChildren(doc.getRootElement());
	}

	/**
	 * For testing purposes
	 * @param args
	 * @throws IOException 
	 * @throws JDOMException 
	 */
	public static void main(String[] args) throws JDOMException, IOException {
		// Auto-generated method stub
		//		XMLInputKit XMLKit = null;
		//		try {
		//			XMLKit = new XMLInputKit();
		//		} catch (IOException | ParserConfigurationException | SAXException e) {
		//			// Auto-generated catch block
		//			e.printStackTrace();
		//		}

		//System.out.print(XMLAuxTools.countChildren(XMLKit.getRootElement()));
		
		// String uniprot = "p08254";
//		String queryURL = "https://www.disgenet.org/api/gda/gene/uniprot/" + uniprot;
		/*
		 * thanks to http://www.xyzws.com/Javafaq/how-to-use-httpurlconnection-post-data-to-web-server/139
		 */
//		String queryParameters = null;
//		try {
//			queryParameters = "source=" + URLEncoder.encode("ALL", "UTF-8") +
//					"&disease_class=" + URLEncoder.encode("C10", "UTF-8") +
//					"&format=" + URLEncoder.encode("xml", "UTF-8");
//		} catch (UnsupportedEncodingException e1) {
//			// Auto-generated catch block
//			e1.printStackTrace();
//		}
		/*
		 * Try this out to solve code 405:
		 * // Attempted FIX20200717CODE405
		 */
//		URL url = null;
//		try {
//			url = new URL(queryURL + "?" + queryParameters);
//		} catch (MalformedURLException e1) {
//			// Auto-generated catch block
//			e1.printStackTrace();
//		}
//		HashMap<String, String> requestMap = new HashMap<String, String>();
//		requestMap.put("Content-Type", "application/xml");

//		HttpURLConnection connection = null;
//		try {
//			connection = HttpConnectionTools.getConnectionFromURL(url);
//		} catch (IOException e1) {
//			// Auto-generated catch block
//			e1.printStackTrace();
//		}
//		/*
//		 * THANKS TO
//		 * duffymo
//		 * AT
//		 * https://stackoverflow.com/questions/1359689/how-to-send-http-request-in-java
//		 * FOR THE ANSWER
//		 */
//
//
//		try {
//			//Create connection
//			// URL url = new URL(queryURL);
//
//			connection = (HttpURLConnection) url.openConnection();
//			connection.setRequestMethod("GET");
//			connection.setRequestProperty("Content-Type", 
//					"application/xml");
//			connection.setUseCaches(false);
//			//connection.setDoOutput(true);
//
//			/*
//			 * test
//			 */
//			//int responseCode = connection.getResponseCode();
//			//System.out.println("Response Code : " + responseCode);
//
//			//Send request
//			//		    DataOutputStream wr = new DataOutputStream (
//			//		        connection.getOutputStream());
//			//		    wr.writeBytes(queryParameters);
//			/*
//			 * Trying out java tutorial's
//			 * approach
//			 */
//			//OutputStreamWriter wr = new OutputStreamWriter (
//			//		connection.getOutputStream());
//			// wr.write(queryParameters); // Attempted FIX20200717CODE405
//			//wr.flush();
//			//wr.close();
//
//			//Get Response  
//			InputStream is = connection.getInputStream();
			//			BufferedReader rd = new BufferedReader(new InputStreamReader(is));
			//			// StringBuilder response = new StringBuilder(); // or StringBuffer if Java version 5+
			//			StringBuffer response = new StringBuffer(); // or StringBuffer if Java version 5+
			//			String line;
			//			while ((line = rd.readLine()) != null) {
			//				response.append(line);
			//				response.append('\r');
			//			}
			//			
			//			System.out.println(response.toString());

			
			//rd.close();
//			connection.disconnect();
//		}
//		catch (Exception e) {
//			// handle exception
//			e.printStackTrace();
//		}
		ArrayList<String> uniprots = ExcelTools.getUniprots("C:UL/miriade/src/MIRIADE_Olink_sorted_biomarkers.xls");
		HashMap<String, Integer> numOfHits = new HashMap<String, Integer>();
		String uniprot;
		int hitCount;
		for (int i = 0; i < uniprots.size(); i++) {
			uniprot = uniprots.get(i);
			hitCount = 0; // Initialize to 0. If there is an exception, 0 will be the default result.
			try {
			hitCount = countUniprotHits(uniprot);
			}
			catch (IOException e) {
				System.out.println("The number of hits for the uniprot " + uniprot + " is: " + 0);
				continue;
			}
			numOfHits.put(uniprot, hitCount);
			System.out.println("The number of hits for the uniprot " + uniprot + " is: " + hitCount);
		}

	}

}
