// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Synthetic example data (R/sfx.R): the presets of synfix(), the columns
// of published tables of Ontario's corrections data (the OTIS tables on
// segregation, restrictive confinement and deaths in custody, a01 to d07,
// and the tables of use of force in correctional institutions), each a
// set of labels, a range of whole numbers or a synthetic identifier
// (syn-00001, ...), and the draws: for each column, n indices drawn as
// R's sample.int(m, n, replace = TRUE) draws them after set.seed(seed)
// (igh::Rng, src/igk.h), so that a seeded table is R's. No record of the
// sources is carried, only their column names and the values a column
// takes.
#include <cstdint>

#include "crimlab.h"
#include "igk.h"
#include "kutil.h"

namespace {

// type 0: the labels of lev, separated by '|'; 1: the whole numbers lo to
// hi; 2: syn-%05d of 1 to max(2, n %/% 4).
struct Col {
  const char *kind, *name;
  int type, lo, hi;
  const char *lev;
};

const Col kCols[] = {
    {"otis-a01", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-a01", "UniqueIndividual_ID", 2, 0, 0,
     ""},
    {"otis-a01", "Region_AtTimeOfPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-a01", "Region_MostRecentPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-a01", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-a01", "Age_Category", 0, 0, 0,
     "18 to 24|25 to 49|50+"},
    {"otis-a01", "MentalHealth_Alert", 0, 0, 0,
     "Yes|No"},
    {"otis-a01", "SuicideRisk_Alert", 0, 0, 0,
     "Yes|No"},
    {"otis-a01", "SuicideWatch_Alert", 0, 0, 0,
     "Yes|No"},
    {"otis-a01", "Number_Of_Placements", 1, 0, 80,
     ""},
    {"otis-b01", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-b01", "UniqueIndividual_ID", 2, 0, 0,
     ""},
    {"otis-b01", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-b01", "Region_AtTimeOfPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-b01", "Region_MostRecentPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-b01", "Age_Category", 0, 0, 0,
     "18 to 24|25 to 49|50+"},
    {"otis-b01", "NumberConsecutiveDays_Segregation", 1, 0, 80,
     ""},
    {"otis-b01", "SegReason_SecurityOfInstitution_SafetyOfOthers", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "SegReason_InmateNeedsProtection", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "SegReason_InmateNeedsProtection_Medical", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "SegReason_SecurityOfInstitution_SafetyOfOthers_Medical", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "SegReason_Disciplinary_Segregation", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "SegReason_InmateRefuseSearch_Scan", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "MentalHealth_Alert", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "SuicideRisk_Alert", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "SuicideWatch_Alert", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "SegReason_Other", 0, 0, 0,
     "Yes|No"},
    {"otis-b01", "Number_Of_Placements", 1, 0, 80,
     ""},
    {"otis-b02", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-b02", "UniqueIndividual_ID", 2, 0, 0,
     ""},
    {"otis-b02", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-b02", "Region_MostRecentPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-b02", "Age_Category", 0, 0, 0,
     "18 to 24|25 to 49|50+"},
    {"otis-b02", "TotalAggregatedDays_Segregation", 1, 0, 80,
     ""},
    {"otis-b03", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-b03", "Region_AtTimeOfPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-b03", "Institution_AtTimeOfPlacement", 0, 0, 0,
     "Algoma Treatment and Remand Centre|Brockville Jail|Central East Correctional Centre|Central Nort"
     "h Correctional Centre|Elgin-Middlesex Detention Centre|Fort Frances Jail|Hamilton Wentworth Dete"
     "ntion Centre|Kenora Jail|Maplehurst Correctional Complex|Monteith Correctional Complex|Niagara D"
     "etention Centre|North Bay Jail|Ontario Correctional Institute|Ottawa-Carleton Detention Centre|Q"
     "uinte Detention Centre|Sarnia Jail|South West Detention Centre|St. Lawrence Valley Correctional "
     "and Treatment Centre|Stratford Jail|Sudbury Jail|Thunder Bay Correctional Centre|Thunder Bay Jai"
     "l|Toronto East Detention Centre|Toronto South Detention Centre|Vanier Centre for Women"},
    {"otis-b03", "Alert_Type", 0, 0, 0,
     "Immigration Hold Flag|Mental Health Alert|Serious Mental Illness Alert|Suicide Risk Alert|Suicid"
     "e Watch Alert|Transgender Alert"},
    {"otis-b03", "Alert_Presence", 0, 0, 0,
     "Yes|No"},
    {"otis-b03", "Number_SegregationPlacements", 1, 0, 80,
     ""},
    {"otis-b04", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-b04", "Region_AtTimeOfPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-b04", "Gender", 0, 0, 0,
     "Female|Male|All"},
    {"otis-b04", "Measure", 0, 0, 0,
     "Maximum|Median|Mode"},
    {"otis-b04", "NumberConsecutiveDays_Segregation", 1, 0, 80,
     ""},
    {"otis-b05", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-b05", "Consecutive_Duration", 0, 0, 0,
     "1 day|2 days|3 days|4 days|5 days|6 to 10 days|11 to 15 days|16 to 20 days|21 to 25 days|26 to 3"
     "0 days|Greater than 30 days"},
    {"otis-b05", "Number_SegregationPlacements", 1, 0, 80,
     ""},
    {"otis-b06", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-b06", "Region_AtTimeOfPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-b06", "Institution_AtTimeOfPlacement", 0, 0, 0,
     "Algoma Treatment and Remand Centre|Brockville Jail|Central East Correctional Centre|Central Nort"
     "h Correctional Centre|Elgin-Middlesex Detention Centre|Fort Frances Jail|Hamilton Wentworth Dete"
     "ntion Centre|Kenora Jail|Maplehurst Correctional Complex|Monteith Correctional Complex|Niagara D"
     "etention Centre|North Bay Jail|Ontario Correctional Institute|Ottawa-Carleton Detention Centre|Q"
     "uinte Detention Centre|Sarnia Jail|South West Detention Centre|St. Lawrence Valley Correctional "
     "and Treatment Centre|Stratford Jail|Sudbury Jail|Thunder Bay Correctional Centre|Thunder Bay Jai"
     "l|Toronto East Detention Centre|Toronto South Detention Centre|Vanier Centre for Women"},
    {"otis-b06", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-b06", "Reason", 0, 0, 0,
     "Close Confinement (Disciplinary Segregation)|Inmate Needs Protection|Inmate Needs Protection: Me"
     "dical|Inmate Refused Search/Scan|Security of Institution/Safety of Others|Security of Institutio"
     "n/Safety of Others: Medical|Other"},
    {"otis-b06", "Number_SegregationPlacements", 1, 0, 80,
     ""},
    {"otis-b07", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-b07", "Alert_Type", 0, 0, 0,
     "Immigration Hold Flag|Mental Health Alert|Serious Mental Illness Alert|Suicide Risk Alert|Suicid"
     "e Watch Alert|Transgender Alert"},
    {"otis-b07", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-b07", "Number_Segregation_Placements_Without_Alert", 1, 0, 80,
     ""},
    {"otis-b07", "Number_Segregation_Placements_With_Alert", 1, 0, 80,
     ""},
    {"otis-b08", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-b08", "Region_AtTimeOfPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-b08", "Institution_AtTimeOfPlacement", 0, 0, 0,
     "Algoma Treatment and Remand Centre|Brockville Jail|Central East Correctional Centre|Central Nort"
     "h Correctional Centre|Elgin-Middlesex Detention Centre|Fort Frances Jail|Hamilton Wentworth Dete"
     "ntion Centre|Kenora Jail|Maplehurst Correctional Complex|Monteith Correctional Complex|Niagara D"
     "etention Centre|North Bay Jail|Ontario Correctional Institute|Ottawa-Carleton Detention Centre|Q"
     "uinte Detention Centre|Sarnia Jail|South West Detention Centre|St. Lawrence Valley Correctional "
     "and Treatment Centre|Stratford Jail|Sudbury Jail|Thunder Bay Correctional Centre|Thunder Bay Jai"
     "l|Toronto East Detention Centre|Toronto South Detention Centre|Vanier Centre for Women"},
    {"otis-b08", "Gender", 0, 0, 0,
     "Female|Male|All"},
    {"otis-b08", "Measure", 0, 0, 0,
     "Median|Mode"},
    {"otis-b08", "NumberConsecutiveDays_Segregation", 1, 0, 80,
     ""},
    {"otis-b09", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-b09", "NumberPlacements_Segregation", 0, 0, 0,
     "1 placement|2 placements|3 placements|4 placements|5 placements|6 to 10 placements|11 to 15 plac"
     "ements|16 to 20 placements|21 to 25 placements|26 to 30 placements|31 to 35 placements|36 to 40 "
     "placements|more than 40 placements"},
    {"otis-b09", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-b09", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c01", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c01", "NumberIndividuals_InCustody", 1, 0, 80,
     ""},
    {"otis-c01", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-c01", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c01", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c02", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c02", "Region_MostRecentPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-c02", "Institution_MostRecentPlacement", 0, 0, 0,
     "Algoma Treatment and Remand Centre|Brockville Jail|Central East Correctional Centre|Central Nort"
     "h Correctional Centre|Elgin-Middlesex Detention Centre|Fort Frances Jail|Hamilton Wentworth Dete"
     "ntion Centre|Kenora Jail|Maplehurst Correctional Complex|Monteith Correctional Complex|Niagara D"
     "etention Centre|North Bay Jail|Ontario Correctional Institute|Ottawa-Carleton Detention Centre|Q"
     "uinte Detention Centre|Sarnia Jail|South West Detention Centre|St. Lawrence Valley Correctional "
     "and Treatment Centre|Stratford Jail|Sudbury Jail|Thunder Bay Correctional Centre|Thunder Bay Jai"
     "l|Toronto East Detention Centre|Toronto South Detention Centre|Vanier Centre for Women"},
    {"otis-c02", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-c02", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c02", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c03", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c03", "Race", 0, 0, 0,
     "Another Race Category|Black|East Asian|Indigenous|Latino|Middle Eastern|More Than One Reported R"
     "ace Category|South Asian|Unknown Or Not Reported|White"},
    {"otis-c03", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-c03", "NumberIndividuals_InCustody", 1, 0, 80,
     ""},
    {"otis-c03", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c03", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c04", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c04", "Region_MostRecentPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-c04", "Race", 0, 0, 0,
     "Another Race Category|Black|East Asian|Indigenous|Latino|Middle Eastern|More Than One Reported R"
     "ace Category|South Asian|Unknown Or Not Reported|White"},
    {"otis-c04", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c04", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c05", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c05", "Region_MostRecentPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-c05", "Religion", 0, 0, 0,
     "Another Religion/Spiritual Affiliation|Buddhist|Christian|Hindu|Indigenous Spirituality|Jewish|M"
     "ore Than One Reported Religion or Spiritual Affiliation|Muslim|No Religion|Sikh|Unknown Or Not R"
     "eported"},
    {"otis-c05", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c05", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c06", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c06", "Region_MostRecentPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-c06", "Age_Category", 0, 0, 0,
     "18 to 24|25 to 49|50+"},
    {"otis-c06", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c06", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c07", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c07", "Alert_Type", 0, 0, 0,
     "Immigration Hold Flag|Mental Health Alert|Serious Mental Illness Alert|Suicide Risk Alert|Suicid"
     "e Watch Alert|Transgender Alert"},
    {"otis-c07", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-c07", "NumberIndividuals_InCustody", 1, 0, 80,
     ""},
    {"otis-c07", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c07", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c08", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c08", "Religion", 0, 0, 0,
     "Another Religion/Spiritual Affiliation|Buddhist|Christian|Hindu|Indigenous Spirituality|Jewish|M"
     "ore Than One Reported Religion or Spiritual Affiliation|Muslim|No Religion|Sikh|Unknown Or Not R"
     "eported"},
    {"otis-c08", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-c08", "NumberIndividuals_InCustody", 1, 0, 80,
     ""},
    {"otis-c08", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c08", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c09", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c09", "Age_Category", 0, 0, 0,
     "18 to 24|25 to 49|50+"},
    {"otis-c09", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-c09", "NumberIndividuals_InCustody", 1, 0, 80,
     ""},
    {"otis-c09", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c09", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c10", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c10", "Region_MostRecentPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-c10", "Institution_MostRecentPlacement", 0, 0, 0,
     "Algoma Treatment and Remand Centre|Brockville Jail|Central East Correctional Centre|Central Nort"
     "h Correctional Centre|Elgin-Middlesex Detention Centre|Fort Frances Jail|Hamilton Wentworth Dete"
     "ntion Centre|Kenora Jail|Maplehurst Correctional Complex|Monteith Correctional Complex|Niagara D"
     "etention Centre|North Bay Jail|Ontario Correctional Institute|Ottawa-Carleton Detention Centre|Q"
     "uinte Detention Centre|Sarnia Jail|South West Detention Centre|St. Lawrence Valley Correctional "
     "and Treatment Centre|Stratford Jail|Sudbury Jail|Thunder Bay Correctional Centre|Thunder Bay Jai"
     "l|Toronto East Detention Centre|Toronto South Detention Centre|Vanier Centre for Women"},
    {"otis-c10", "Gender", 0, 0, 0,
     "Female|Male|All"},
    {"otis-c10", "Measure", 0, 0, 0,
     "Maximum|Median|Mode"},
    {"otis-c10", "TotalAggregatedDays_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c10", "TotalAggregatedDays_Segregation", 1, 0, 80,
     ""},
    {"otis-c11", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c11", "Aggregate_Duration", 0, 0, 0,
     "1 day|2 days|3 days|4 days|5 days|6 to 10 days|11 to 15 days|16 to 20 days|21 to 25 days|26 to 3"
     "0 days|Greater than 30 days"},
    {"otis-c11", "NumberIndividuals_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c11", "NumberIndividuals_Segregation", 1, 0, 80,
     ""},
    {"otis-c12", "EndFiscalYear", 1, 2018, 2024,
     ""},
    {"otis-c12", "Region_MostRecentPlacement", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-c12", "Gender", 0, 0, 0,
     "Female|Male|All"},
    {"otis-c12", "Measure", 0, 0, 0,
     "Maximum|Median|Mode"},
    {"otis-c12", "TotalAggregatedDays_RestrictiveConfinement", 1, 0, 80,
     ""},
    {"otis-c12", "TotalAggregatedDays_Segregation", 1, 0, 80,
     ""},
    {"otis-d01", "Year", 1, 2018, 2024,
     ""},
    {"otis-d01", "UniqueIndividual_ID", 2, 0, 0,
     ""},
    {"otis-d01", "Region_AtTimeOfDeath", 0, 0, 0,
     "Central|Eastern|Northern|Toronto|Western"},
    {"otis-d01", "HousingUnit_Type", 0, 0, 0,
     "General Population|Segregation|Specialized Care|Protective Custody Outside of a Facility"},
    {"otis-d01", "MedicalCauseOfDeath", 0, 0, 0,
     "Drug Toxicity|Natural Causes|Other|Unknown / To Be Determined"},
    {"otis-d01", "MeansOfDeath", 0, 0, 0,
     "Accident|Homicide|Natural Causes|Suicide|Undetermined|To Be Determined"},
    {"otis-d02", "Year", 1, 2018, 2024,
     ""},
    {"otis-d02", "Gender", 0, 0, 0,
     "Female|Male"},
    {"otis-d02", "Number_CustodialDeaths", 1, 0, 80,
     ""},
    {"otis-d03", "Year", 1, 2018, 2024,
     ""},
    {"otis-d03", "Race", 0, 0, 0,
     "Another Race Category Black|East Asian|Indigenous|Latino|Middle Eastern|More Than One Reported R"
     "ace Category|South Asian|Unknown Or Not Reported|White"},
    {"otis-d03", "Number_CustodialDeaths", 1, 0, 80,
     ""},
    {"otis-d04", "Year", 1, 2018, 2024,
     ""},
    {"otis-d04", "Religion", 0, 0, 0,
     "Another Religion/Spiritual Affiliation|Buddhist|Christian|Hindu|Indigenous Spirituality|Jewish|M"
     "uslim|More Than One Reported Religion or Spiritual Affiliation|No Religion|Sikh|Unknown Or Not R"
     "eported"},
    {"otis-d04", "Number_CustodialDeaths", 1, 0, 80,
     ""},
    {"otis-d05", "Year", 1, 2018, 2024,
     ""},
    {"otis-d05", "Age_Category", 0, 0, 0,
     "18 to 24|25 to 49|50+"},
    {"otis-d05", "Number_CustodialDeaths", 1, 0, 80,
     ""},
    {"otis-d06", "Year", 1, 2018, 2024,
     ""},
    {"otis-d06", "Alert_Type", 0, 0, 0,
     "Mental Health Alert|Suicide Risk Alert|Suicide Watch Alert"},
    {"otis-d06", "MedicalCauseOfDeath", 0, 0, 0,
     "Drug Toxicity|Natural Causes|Other|Unknown / To Be Determined"},
    {"otis-d06", "Number_CustodialDeaths", 1, 0, 80,
     ""},
    {"otis-d07", "Year", 1, 2018, 2024,
     ""},
    {"otis-d07", "Alert_Type", 0, 0, 0,
     "Mental Health Alert|Suicide Risk Alert|Suicide Watch Alert"},
    {"otis-d07", "HousingUnit_Type", 0, 0, 0,
     "General Population / Protective Custody Unit|Outside of a Facility|Segregation Conditions|Specia"
     "lized Care Unit"},
    {"otis-d07", "Number_CustodialDeaths", 1, 0, 80,
     ""},
    {"cuof-incidents", "year", 1, 2023, 2024,
     ""},
    {"cuof-incidents", "incident_ID", 2, 0, 0,
     ""},
    {"cuof-incidents", "status", 0, 0, 0,
     "Final"},
    {"cuof-incident_type", "year", 1, 2023, 2024,
     ""},
    {"cuof-incident_type", "incident_ID", 2, 0, 0,
     ""},
    {"cuof-incident_type", "incident_type", 0, 0, 0,
     "Use of Force|Assault|Threat"},
    {"cuof-indigenous", "year", 1, 2023, 2024,
     ""},
    {"cuof-indigenous", "inmate_ID", 2, 0, 0,
     ""},
    {"cuof-indigenous", "force_used", 0, 0, 0,
     "Yes|Unknown"},
    {"cuof-inmate_incident", "year", 1, 2023, 2024,
     ""},
    {"cuof-inmate_incident", "incident_ID", 2, 0, 0,
     ""},
    {"cuof-inmate_incident", "inmate_ID", 2, 0, 0,
     ""},
    {"cuof-inmate_participant", "year", 1, 2023, 2024,
     ""},
    {"cuof-inmate_participant", "inmate_ID", 2, 0, 0,
     ""},
    {"cuof-inmate_participant", "force_used", 0, 0, 0,
     "Yes|Unknown"},
    {"cuof-institution_summary", "year", 1, 2023, 2024,
     ""},
    {"cuof-institution_summary", "region", 0, 0, 0,
     "Central|Eastern"},
    {"cuof-institution_summary", "institution", 0, 0, 0,
     "Hamilton-Wentworth Detention Centre|Maplehurst Correctional Complex|Niagara Detention Centre|Van"
     "ier Centre for Women|Brockville Jail"},
    {"cuof-institution_summary", "number_of_incidents", 1, 36, 624,
     ""},
    {"cuof-institution_summary", "number_unique_inmate_participant", 1, 35, 758,
     ""},
    {"cuof-institution_summary", "number_unique_staff_participant", 1, 62, 429,
     ""},
    {"cuof-institution_summary", "use_of_force_rate", 1, 11, 24,
     ""},
    {"cuof-location_summary", "year", 1, 2023, 2024,
     ""},
    {"cuof-location_summary", "general_location", 0, 0, 0,
     "Administration Areas/Offices|Courthouse|Elevator|Hallway/Stairwell|Health Care Area"},
    {"cuof-location_summary", "number_of_incidents", 1, 1, 177,
     ""},
    {"cuof-race", "year", 1, 2023, 2024,
     ""},
    {"cuof-race", "inmate_ID", 2, 0, 0,
     ""},
    {"cuof-race", "force_used", 0, 0, 0,
     "Yes|Unknown"},
    {"cuof-ethnic_origin", "year", 1, 2023, 2024,
     ""},
    {"cuof-ethnic_origin", "inmate_ID", 2, 0, 0,
     ""},
    {"cuof-ethnic_origin", "force_used", 0, 0, 0,
     "Yes|Unknown"},
    {"cuof-religion", "year", 1, 2023, 2024,
     ""},
    {"cuof-religion", "inmate_ID", 2, 0, 0,
     ""},
    {"cuof-religion", "force_used", 0, 0, 0,
     "Yes|Unknown"},
    {"cuof-select_incident_summary", "year", 1, 2023, 2024,
     ""},
    {"cuof-select_incident_summary", "incident_type", 0, 0, 0,
     "Assault"},
    {"cuof-select_incident_summary", "level_one", 0, 0, 0,
     "(P1) Serious Inmate on Inmate|(P1) Serious Inmate on Staff (Requires contact)|(P2) Inmate on Inm"
     "ate - No serious injury or hospitalization|(P2) Inmate on Member of Public|(P2) Inmate on Staff "
     "- No Serious Injury or Hospitalization - No contact"},
    {"cuof-select_incident_summary", "number_of_incidents", 1, 3, 1041,
     ""},
    {"cuof-staff_incident", "year", 1, 2023, 2024,
     ""},
    {"cuof-staff_incident", "incident_ID", 2, 0, 0,
     ""},
    {"cuof-staff_incident", "staff_ID", 2, 0, 0,
     ""},
};

}  // namespace

extern "C" int clab_sfxcol(int j, const char **kind, const char **name, int *type, int *lo, int *hi,
                           const char **lev) {
  if (j < 0 || j >= static_cast<int>(sizeof kCols / sizeof kCols[0])) return CLAB_EARG;
  const Col &c = kCols[j];
  *kind = c.kind;
  *name = c.name;
  *type = c.type;
  *lo = c.lo;
  *hi = c.hi;
  *lev = c.lev;
  return CLAB_OK;
}

extern "C" int clab_synfix(const double *dn, int k, int64_t n, uint32_t seed, int *code) CLAB_TRY {
  if (k < 1 || n < 1) return CLAB_EARG;
  for (int j = 0; j < k; ++j) {
    if (!(dn[j] >= 1 && dn[j] <= 2147483647.0)) return CLAB_EARG;
  }
  igh::Rng g(seed);
  for (int j = 0; j < k; ++j) {
    int *c = code + static_cast<int64_t>(j) * n;
    if (kutil::stop(static_cast<double>(n))) return CLAB_EINTR;
    for (int64_t i = 0; i < n; ++i) c[i] = static_cast<int>(g.index(dn[j]));
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
